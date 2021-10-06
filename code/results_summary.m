clear
% This script generates the numerical robustness presented in the paper
% The adjustable parameters may be used the improve our inner and outer aproximations
% In order to obtain the values for a particular process, uncomment the corresponding line

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%  Adjustable parameters  %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n = 13; %This corresponds to the n of the family of uniform qubit states with N=2n^2 states
states = uniform_qubit_states(n); % Pick N=2*n^2 rather uniform qubit states as in Ref. 
% states=rNRandomPureStates(2,2*n^2); %Pick N=2*n^2 uniformly random qubit states as in 
k = 1;  %Setting k=1 is equivalent of PPT, k>=2 also imposes k-sym Bose extention

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%  Adjustable parameters  %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

non_ccdc_processes_generator; % This command generates all process mentioned in the paper.
%Uncomment the line corresponding to the desired process

 ccdc_robustness_summary(W222, [2 2 2], states, k);   
% ccdc_robustness_summary(W224, [2 2 4], states, k); 
% ccdc_robustness_summary(W_MRSR, [2 2 2],states, k);  
% ccdc_robustness_summary(W_PPT, [2 2 2], states, k);   
% ccdc_robustness_summary(W_SEP, [2 2 2], states, k)
%ccdc_robustness_summary(W_FB, [2 2 2 2 8], states); %Be careful, this may take a long time


%states = NRandomPureStates(3,2*n^2); %Since W339 next ones has qutrits in Alice's input, the variable 'states' should be qutrits
%ccdc_robustness_summary(W339, [3 3 9], states, k);  
