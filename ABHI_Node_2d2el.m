classdef ABHI_Node_2d2el < RC_Node_2d1el
% Replace XYZ by your initials and rename the file accordingly before proceeding

% Node class for 2nd order analysis of a 2-dimensional structure
    
    % Protected properties go here
    properties (Access = protected)
        
    end
    
    % Public methods go here
    methods (Access = public)
        function self=ABHI_Node_2d2el(node_num,node_coord)
            self=self@RC_Node_2d1el(node_num,node_coord);
        end
        %% Update Node Coord
        % update the coordinates of the node object after each iteration
        % step by an amount dDelta. dDelta is a 1x2 vector that represents
        % the change in displacement for that step
        
        function UpdateNodeCoord(self,dDelta)
            self.node_coord = self.node_coord + dDelta';    
        end
    end
    
    
    % Protected methods go here
    methods (Access = protected)
        
    end
end
