%% %%%%%%%%%%%%%%%%%%%%%%% RUN ANALYZEPRF

% face semi-automated model using the full data
for p=1:1, nsdfaceprf(p,1,2,1,1:3);, end

% contrastgrid model using the full data
for p=1:8, nsdfaceprf(p,3,1,1,1:3);, end
for p=1:8, nsdfaceprf(p,4,1,1,1:3);, end  % NEW

% face automated model using the full data
for p=1:1, nsdfaceprf(p,1,1,1,1:3);, end

% word semi-automated model using the full data
for p=1:1, nsdfaceprf(p,2,2,1,1:3);, end

% word automated model using the full data
for p=1:1, nsdfaceprf(p,2,1,1,1:3);, end

% at this point ran all 8 subjects maybe?

% proceed to 1.8mm
for p=1:8, nsdfaceprf(p,3,1,1,4);, end
for p=1:8, nsdfaceprf(p,4,1,1,4);, end

% split half!!! face auto
for p=1:8, nsdfaceprf(p,1,1,2,1:3);, end
for p=1:8, nsdfaceprf(p,1,1,3,1:3);, end

% bodies model (automated)
for p=1:8, nsdfaceprf(p,5,1,1,1:3);, end

% foreground model (automated)
for p=1:8, nsdfaceprf(p,6,1,1,1:3);, end
for p=1:8, nsdfaceprf(p,7,1,1,1:3);, end
for p=1:8, nsdfaceprf(p,8,1,1,1:3);, end

% foreground model (automated) 1.8mm volume version both positive and negative
for p=1:8, nsdfaceprf(p,6,1,1,4);, end
for p=1:8, nsdfaceprf(p,7,1,1,4);, end

% contrastgrid model fixed exponent and baseline
for p=1:8, nsdfaceprf(p,9,1,1,1:3);, end

% contrastgrid model using contrastgridNEW preparation
for p=1:8, nsdfaceprf(p,10,1,1,1:3);, end
for p=1:8, nsdfaceprf(p,11,1,1,1:3);, end

% face automated NEGATIVE  [but these results are weird because of offset issues?]
for p=1:8, nsdfaceprf(p,12,1,1,1:3);, end

% cerebellum  [did tk find these interesting?]
for p=1:8
  nsdfaceprf(p,1,1,1,5);   % face auto
  nsdfaceprf(p,2,1,1,5);   % word auto
  nsdfaceprf(p,5,1,1,5);   % bodies
  nsdfaceprf(p,6,1,1,5);   % foreground
  nsdfaceprf(p,8,1,1,5);   % background
  nsdfaceprf(p,10,1,1,5);  % contrastgridNEW
  nsdfaceprf(p,13,1,1,5);  % salience $$$
end

% salience
for p=1:8, nsdfaceprf(p,13,1,1,1:3);, end

% salience, 1.8mm volume version both positive and negative  $$$
for p=1:8, nsdfaceprf(p,13,1,1,4);, end
for p=1:8, nsdfaceprf(p,14,1,1,4);, end

% facebody model (automated)
for p=1:8, nsdfaceprf(p,15,1,1,1:3);, end

%% %%%%%%%%%%%%%%%%%%%%%%% POST-PROCESS THE RESULTS (this includes saving various files)

% category model (note that rr is [] since we save all rrs together. also, this is a hack to help maintain old results)
for p=1:8
  nsdfaceprffinish(p,1,1,[],1:3,1);
  nsdfaceprffinish(p,1,2,[],1:3,1);
  nsdfaceprffinish(p,2,1,[],1:3,1);
  nsdfaceprffinish(p,2,2,[],1:3,1);
end

% fullprf model
for p=1:1
%  nsdfaceprffinish(1,1,2,1,1:3);
%  nsdfaceprffinish(1,3,1,1,1:3);
%  nsdfaceprffinish(1,1,1,1,1:3);
  nsdfaceprffinish(1,2,1,1,1:3);  
  nsdfaceprffinish(1,2,2,1,1:3);
end
for p=2:8
%  nsdfaceprffinish(p,1,1,1,1:3);
%  nsdfaceprffinish(p,2,1,1,1:3);
%  nsdfaceprffinish(p,1,2,1,1:3);
%  nsdfaceprffinish(p,2,2,1,1:3);
  nsdfaceprffinish(p,3,1,1,1:3);
end
for p=1:8
  nsdfaceprffinish(p,4,1,1,1:3);
end
for p=1:8
  nsdfaceprffinish(p,3,1,1,4);
  nsdfaceprffinish(p,4,1,1,4);
end
for p=1:8
  nsdfaceprffinish(p,1,1,2,1:3,[],0);
  nsdfaceprffinish(p,1,1,3,1:3);
end

% NEXT SET: ran on July 1 2020.
for p=1:8
  nsdfaceprffinish(p,5,1,1,1:3);
end
for p=1:8   % reran this to reinstute standard method
  nsdfaceprffinish(p,3,1,1,1:3);
end

% NEXT SET: ran on July 5+12 2020.
for p=1:8
  nsdfaceprffinish(p,6,1,1,1:3);
  nsdfaceprffinish(p,8,1,1,1:3);
end

% NEXT SET: ran on July 9 2020.
for p=1:8
  nsdfaceprffinish(p,6,1,1,4);
  nsdfaceprffinish(p,7,1,1,4);
end

% july 21 2020:
for p=1:8
  nsdfaceprffinish(p,7,1,1,1:3);
end

% sep 19 2020:
for p=1:8
  nsdfaceprffinish(p,9,1,1,1:3);
end

% sep 23-26 2020:
for p=1:8
  nsdfaceprffinish(p,10,1,1,1:3);
  nsdfaceprffinish(p,11,1,1,1:3);
end

% oct 4 2020:
for p=1:8
  nsdfaceprffinish(p,12,1,1,1:3);
end

% oct 13 2020
for p=1:8
  for wh=[1 2 5 6 8 10]
    nsdfaceprffinish(p,wh,1,1,5);
  end
end

% nov 1 2020:
for p=1:8
  nsdfaceprffinish(p,13,1,1,1:3);
end

% jun 30 2021:
for p=1:8
  nsdfaceprffinish(p,15,1,1,1:3);
end

% record of what we have:
% categoryresults_subj01-8.mat   - simple face vs non-face
% XXX prfresults_subj01-8.mat    - quick mode fitting for faces and words
% fullprfresults_subj01-8.mat    - full pRF fit for contrast grid model, face/word models
