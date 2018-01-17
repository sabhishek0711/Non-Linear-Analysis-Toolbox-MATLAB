classdef ABHI_Element_2d2el < RC_Element_2d1el
% Replace XYZ by your initials and rename the file accordingly before proceeding

% Element class for 2nd order analysis of a 2-dimensional structure
    
    % Protected properties go here
    properties (Access = protected)
       k_local
       kg_local
       
       % global force vector of the element used for computing error
       f_global
       
    end
    
    % Public methods go here
    methods (Access = public)
        function self=ABHI_Element_2d2el(element_nodes, A, Ayy, Izz, E, v, truss)
            self=self@RC_Element_2d1el(element_nodes, A, Ayy, Izz, E, v, truss);
            self.f_local = zeros(6,1); %set local force vector to zero to start
            self.ComputeElementGeomStiffMatrix();
        end
        
       %% Compute Forces 
        function ComputeForces(self, del_global)
            self.del_global = del_global;
            
            % Compute the element displacement vector in local coordinates
            self.del_local = self.gamma * self.del_global;
            
            % Compute the element force vector in local coordinates
            %self.f_local = self.ke_local * self.del_local;
            
            un=(self.del_local(4)-self.del_local(1)) + ((self.del_local(5)-self.del_local(2))^2 + (self.del_local(4)-self.del_local(1))^2)/(2*self.L);
            thetaR=atan((self.del_local(5)-self.del_local(2))/(self.L + self.del_local(4)-self.del_local(1)));
            thetaAN=self.del_local(3)-thetaR;
            thetaBN=self.del_local(6)-thetaR;
            deltaN=[0;0;thetaAN;un;0;thetaBN];
            dF = (self.ke_local + self.kg_local)*deltaN;
            self.f_local= self.f_local + dF;
        end
        %% Compute Global Stiffness        
        function ComputeGlobalStiffnessMatrix(self)
            self.k_local=self.ke_local + self.kg_local;
            self.k_global=self.gamma'*self.k_local*self.gamma;  
        end
        
        %% Update Transformation Matrix
        % This function updates the transformation matrix of the element
        function UpdateTransformationMatrix(self)
            self.ComputeTransformationMatrix()
        end
        %% Compute Element Geometric Stiffness Matrix
        % Computes the geoemetric stiffness matrix for the element
        function ComputeElementGeomStiffMatrix(self)
            
            %self.ComputeGlobalStiffnessMatrix();
            %self.Ke_global=self.k_global;
            %Ke_Local=self.ke_local;
            L=self.L;
            P = self.f_local(4);
            
            
            Kg=(P/L)*[1  0   0  -1   0   0;...
                      0 6/5 L/10 0 -6/5 L/10;...
                      0 L/10 2*L^2/15 0 -L/10 -L^2/30;...
                      -1  0   0   1   0   0;...
                      0 -6/5 -L/10 0 6/5 -L/10;...
                      0 L/10 -L^2/30 0 -L/10 2*L^2/15];
            
            self.kg_local=sparse(Kg);      
        end
        
        %% Compute Global Element Forces
        % computes the elements global forces in global coordinates based
        % on updated geometry for purposes of error calculation
        function ComputeGlobalElementForces(self)
            self.f_global = self.gamma'*self.f_local;
        end
        
        %% Get F global
        % retrieves the elements global internal forces for access in the
        % analysis class error calculation
        function f_global = GetFGlobal(self)
            f_global = self.f_global;
        end    
    end
    
    % Protected methods go here
    methods (Access = protected)
           
    end
end
