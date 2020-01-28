% This script iterates witness sdp for a given initial process matrix W_0
% with causalmaxviol, in order to find a process that violates witness
% significantly
% This code works for either order problem or ccdc problem

function [J,Jf,W] = noncausal_processfinder(d,W_0,robustness,probtype)

tic

W = cell(1);
J = [];
i=1;

if isequal(probtype,'order')
    while 1
    [Ji,Si] = order_sdp(d,robustness,'primal',W_0);
    J = [J Ji];
 
    [Wi,Jk] = causalwitmaxviol(d,Si);
    W_0 = Wi;
    W{i} = W_0;
    J = [J Jk];
    
    i=i+1
    Jf = J(end)
    
    if isequal(robustness,'GR')
        if Jf < -0.2
            break
        end  
    else
        if Jf < -1.5
            break
        end
    end
           
        
    
    if abs(J(end) - J(end-1)) < 0.00001
       break 
    end
        
    t = toc;
    disp(datestr(datenum(0,0,0,0,0,t),'HH:MM:SS'))        
    end
    
elseif isequal(probtype,'ccdc')

    while 1
    [Ji,Si] = ccdc_complete(d,robustness,'primal',W_0);
    J = [J Ji];
 
    [Wi,Jk] = ccdcwitmaxviol(d,Si);
    W_0 = Wi;
    W{i} = W_0;
    J = [J Jk];
    
    i=i+1
    Jf = J(end)
    
%     if isequal(robustness,'GR')
%         if Jf < -0.2
%             break
%         end  
%     else
%         if Jf < -1.5
%             break
%         end
%     end
           
        
    
    if abs(J(end) - J(end-1)) < 0.0000001
       break 
    end
        
    t = toc;
    disp(datestr(datenum(0,0,0,0,0,t),'HH:MM:SS'))        
    end

end