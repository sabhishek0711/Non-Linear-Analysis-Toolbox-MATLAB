classdef ABHI_Analysis_2d2el < RC_Analysis_2d1el
% Replace XYZ by your initials and rename the file accordingly before proceeding

% Analysis class for 2nd order analysis of a 2-dimensional structure
    
    % Protected properties go here
    properties (Access = protected)
        numsteps
        
        % maxsteps is the maximum number of steps that code will run
        % through. It is determined by the lesser of numsteps or
        % stop_ratio/ratio_req
        maxsteps
        
        ratio_req
        stop_ratio
        restart
        apratios
        limit_state
        h_stat_mes
        
        %dconcen and dconcen_t are the incremental concentrated loads
        %stored in addition to the full concentrated loads
        dconcen
        dconcen_t
        
        %step is an integer indicating what the current step of the
        %analysis is
        step
        
        % R is a vector that is ndofx1 to store the calculated reactions after each step
        R
        
        % E is a matrix of dimenstions free_dof x m that stores the computed error
        % vector for the mth step
        E
        
        %vectors that are mx1 to store the energy indices after for the mth
        %step
        load_norm_E
        energy_norm_E
        
        % DEFL_STEP and REACT_STEP are the equivalent of
        % DEFL and REACT except they only store values for the
        % current iteration and are updated every iteration
        DEFL_STEP
        REACT_STEP
        
        % New output values see ud_2d2el
        APRATIOS
        LIMIT_STATE
    end

    % Public methods go here
    methods (Access = public)
        function self=ABHI_Analysis_2d2el(nnodes, coord, fixity, concen,...
                nele, ends, A, Ayy, Izz, E, v, truss,numsteps,ratio_req,...
                stop_ratio,restart,apratios,limit_state,h_stat_mes)
            
            self=self@RC_Analysis_2d1el(nnodes, coord, fixity, concen, nele, ends, A, Ayy, Izz, E, v, truss);
            
            self.numsteps=numsteps;
            self.maxsteps=min(numsteps,floor(stop_ratio/ratio_req));
            self.ratio_req=ratio_req;
            self.stop_ratio=stop_ratio;
            self.restart=restart;
            self.apratios=apratios;
            self.limit_state=limit_state;
            self.h_stat_mes=h_stat_mes;
            self.dconcen = concen.*ratio_req;
            self.dconcen_t = self.dconcen';
            
            
        end
        %% Run Analysis
        % Run 2nd order analysis
        function RunAnalysis(self)
            self.InitializeOutputVariables()
            self.CreateStiffnessMatrix();
            self.CreateLoadVectors();
            self.step = 1;
            
            
            % Run the analysis only if the structure is stable, i.e. the Kff matrix is well conditioned
            if self.AFLAG
                while self.step <= self.maxsteps && self.LIMIT_STATE == 0
                   % Append the applied load ratio for the current step to
                   % the APRATIOS vector
                    self.APRATIOS = [self.APRATIOS;self.step*self.ratio_req];
                    self.runIteration()
                    self.step = self.step + 1;
                end
            end
            
            %Plot error outpus
            RC_Plot_Errors(self.load_norm_E,self.energy_norm_E,self.APRATIOS);
            
        end
        
        %% Run Iteration
        % runs each iteration of the 2nd order analysis
        function runIteration(self)
%             text_mess = ['Performing step #',num2str(self.step)];
%             set(h_stat_mes,'String',text_mess); drawnow;
            self.ComputeDisplacementsReactions(self.step);
            self.RecoverElementForces(self.step) 
            self.ChangeNodes();
            self.UpdateStiffnessMatrices();
            self.CreateStiffnessMatrix();
            self.CheckLimitState();
            self.ComputeError();
        end
        
        %% Get Mastan2 Returns
        %  Returns the matrices that need to be returned to Mastan2
        function [DEFL, REACT, ELE_FOR, AFLAG, APRATIOS, LIMIT_STATE] = GetMastan2Returns(self)
            DEFL = self.DEFL;
            REACT = self.REACT;
            ELE_FOR = self.ELE_FOR;
            AFLAG = self.AFLAG;
            APRATIOS = self.APRATIOS;
            LIMIT_STATE = self.LIMIT_STATE; 
        end
       
    end
    
    % Protected methods go here
    methods (Access = protected)
        %% Create Nodes
         %Create the nnodes x 1 vector of node objects representing all
         %the nodes in the structure as second order object
        function CreateNodes(self)
            for i = 1:self.nnodes
                
                % Create a Node object and append it to the "nodes" vector
                self.nodes = [self.nodes; ABHI_Node_2d2el(i, self.coord_t(:,i))];
            end
        end
        
        %% Create Elements
        %  Create the nele x 1 vector of element objects representing all the elements in the structure
        %    Arguments are all matrices received from Mastan2. Refer to comments in ud_2d1el.m for details.
        function CreateElements(self, A, Ayy, Izz, E, v)
            for i = 1:self.nele
                % Create an Element object and append it to the "elements" vector
                self.elements = [self.elements; ABHI_Element_2d2el(self.nodes(self.ends(i, 1:2)), A(i), ...
                                    Ayy(i), Izz(i), E(i), v(i), self.truss)];
            end
        end
        %% Initialize Output Variables
        function InitializeOutputVariables(self)
            % Note: DEFL_STEP is 3xnnodes rather than nnodesx3 like the
            % output variable DEFL. This comes into play when using
            % DEFL_STEP for linear indexing. REACT_STEP is the same but
            
            self.DEFL_STEP = zeros(self.num_dof_node, self.nnodes);
            self.REACT_STEP = zeros(self.num_dof_node, self.nnodes);
            
            
            self.DEFL = zeros(self.nnodes,self.num_dof_node);
            self.REACT = zeros(self.nnodes,self.num_dof_node);
            self.ELE_FOR = zeros(self.nele,self.num_dof_node*2);

            self.APRATIOS = self.apratios; %initialize applied load ratios
            self.E = []; %initialize the error vector
            self.load_norm_E = [];
            self.energy_norm_E = [];
            self.LIMIT_STATE = 0; %set equal to zero for first iteration
           
        end
        
        %% Create Load Vectors
        %  Create the applied load vectors
        function CreateLoadVectors(self)
            
            % Compute vector of concentrated loads applied at the free and support degrees of freedom using
            % linear indexing of the "concen_t" matrix
            self.Pf = self.dconcen_t(self.dof_free);
            self.Psupp = self.dconcen_t(self.dof_supp);
            
            % Compute vector of specified displacements using linear indexing of the "fixity_t" matrix
            self.deln = self.fixity_t(self.dof_disp);
            
        end
        
        %% Compute Displacements Reactions
        %  Compute the displacements and reactions and format them to return to Mastan2
        function ComputeDisplacementsReactions(self,i)
            
            % Compute the displacements
            self.delf = self.Kff \ (self.Pf - self.Kfn*self.deln);
            
            % Compute the reactions, accounting for loads applied directly on the supports
            self.Ps = self.Ksf*self.delf + self.Ksn*self.deln - self.Psupp;
            self.Pn = self.Knf*self.delf + self.Knn*self.deln;
            
            % Format the computed displacements using linear indexing of the "DEFL" matrix
            self.DEFL_STEP(self.dof_free) = self.delf;
            self.DEFL_STEP(self.dof_disp) = self.deln;
            if i == 1
                self.DEFL(:,:,i) = self.DEFL_STEP';
            else
                self.DEFL = cat(3,self.DEFL,self.DEFL(:,:,i-1)+self.DEFL_STEP');
            end
            
            
            % Format the computed reactions using linear indexing of the "REACT" matrix
            self.REACT_STEP(self.dof_supp) = self.Ps;
            self.REACT_STEP(self.dof_disp) = self.Pn;
            if i == 1
                self.REACT(:,:,i) = self.REACT_STEP';
            else
                %self.REACT(:,:,i) = self.REACT(:,:,i-1)+self.REACT_STEP';
                self.REACT = cat(3,self.REACT,self.REACT(:,:,i-1)+self.REACT_STEP');
            end    
        end
        
        %% Recover Element Forces
        %  Recover the local element forces and format them to return to Mastan2
        function RecoverElementForces(self, j)
            % add another layer filled with zeros to the ELE_FOR matrix if not first iteration
            if j ~= 1
            self.ELE_FOR = cat(3,self.ELE_FOR,zeros(self.nele,self.num_dof_node*2));
            end
            for i = 1:self.nele
                
                % Obtain the displacements at the degrees of freedom corresponding to element i using linear
                % indexing of the "DEFL_STEP" matrix
                if j==1
                self.elements(i).ComputeForces(self.DEFL_STEP(self.elements(i).GetElementDOF()));
                else
                self.elements(i).ComputeForces(self.DEFL_STEP(self.elements(i).GetElementDOF()));    
                end
                %add in calculated F Local forces for element i to the
                %ELE_FOR matrix. i is element number j is the step number.
                self.ELE_FOR(i,:,j) = self.elements(i).GetFLocal();
            end
        end
        
        %% Change Nodes
        % change the coordinates of the nodes in the structure based off of
        % the deflections calculated in the previous step
        function ChangeNodes(self)
            %loop through all nodes to change their coordinates
            for i = 1:self.nnodes
                DEFL_STEP_t = self.DEFL_STEP';
                %Get the deflections at each node for that step
                self.nodes(i).UpdateNodeCoord(DEFL_STEP_t(i,[1,2]))
            end
        end
        

        %% Update Stiffness Matrices
        % update the gamma and kg matrices for the next iteration
        function UpdateStiffnessMatrices(self)
            for i = 1:self.nele
                self.elements(i).UpdateTransformationMatrix();
                self.elements(i).ComputeElementGeomStiffMatrix();
            end
        end
        
        %% Compute Error
        % Computes error vector for all free degrees of freedom
        function ComputeError(self)
            %R is a ndofx1 vector to store the summed element forces at
            %each degree of freedom. The reactions at the free dof are
            %extracted at the end
            self.R = zeros(self.num_dof_total,1); 
            for i = 1:self.nele
                self.elements(i).ComputeGlobalElementForces(); 
                f_global = self.elements(i).GetFGlobal;
                % vector containing the number of the free degrees of
                % freedom for that element
                element_dof = self.elements(i).GetElementDOF();
               
               % assemble R vector every degree of freedom
               self.R(element_dof) = self.R(element_dof) + f_global;
            end
            applied_load = self.APRATIOS(self.step)*self.concen_t(:);
            E_total = applied_load-self.R;
            
            %The error vector for the current step at free dof
            E_step = E_total(self.dof_free); 
            
            % append to the E vector from the previous step. 
            self.E = horzcat(self.E,E_step);  
            
            % compute the error indices for the current step
            
            % Load norm error index
            applied_load_free = applied_load(self.dof_free);
            load_norm_E = norm(E_step)/norm(applied_load_free);
            self.load_norm_E = [self.load_norm_E;load_norm_E];
            
            % Energy norm error index
            energy_norm_E = (abs(E_step)'*abs(self.delf))/(abs(applied_load_free)'*abs(self.delf));
            self.energy_norm_E = [self.energy_norm_E;energy_norm_E];
        end
            
        %% Check Limit State
        % checks what the limit state of the current Kff matrix is
        function CheckLimitState(self)
            [~, p] = chol(self.Kff);
            if p ==0
                %structure is loading
                self.LIMIT_STATE = 0;
            else
                %structure is unloading
                self.LIMIT_STATE = 1;
            end
        end
    end
end

