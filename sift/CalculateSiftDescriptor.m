
function [database, lenStat] = CalculateSiftDescriptor(rt_img_dir, rt_data_dir, gridSpacing, patchSize, maxImSize, nrml_threshold)
%==========================================================================
% usage: calculate the sift descriptors given the image directory
%
% inputs
% rt_img_dir    -image database root path
% rt_data_dir   -feature database root path
% gridSpacing   -spacing for sampling dense descriptors
% patchSize     -patch size for extracting sift feature
% maxImSize     -maximum size of the input image
% nrml_threshold    -low contrast normalization threshold
%
% outputs
% database      -directory for the calculated sift features
%
% Lazebnik's SIFT code is used.
%
% written by Jianchao Yang
% Mar. 2009, IFP, UIUC
%==========================================================================


disp('Extracting SIFT features...');
subfolders = dir(rt_img_dir);

siftLens = [];

database = [];

database.imnum = 0; % total image number of the database
database.cname = {}; % name of each class
database.label = []; % label of each class
database.path = {}; % contain the pathes for each image of each class
database.nclass = 0;

for ii = 1:length(subfolders)
    subname = subfolders(ii).name;
    
    if ~strcmp(subname, '.') & ~strcmp(subname, '..')
        database.nclass = database.nclass + 1;
        
        database.cname{database.nclass} = subname;
        
        frames = dir(fullfile(rt_img_dir, subname, '*.jpg'));
        
        c_num = length(frames);           
        database.imnum = database.imnum + c_num;
        database.label = [database.label; ones(c_num, 1)*database.nclass];
        
        siftpath = fullfile(rt_data_dir, subname);        
        if ~isdir(siftpath)
            mkdir(siftpath);
        end;
        
        for jj = 1:c_num
            imgpath = fullfile(rt_img_dir, subname, frames(jj).name);
            
            I = imread(imgpath);
            if ndims(I) == 3
                I = single(im2double(rgb2gray(I)));
            else
                I = single(im2double(I));
        
            end;
           if ndims(I) == 3
%                  I = rgb2gray(I);
% I = single(rgb2gray(I)) ;
            else
                I = I;
            end;
%              I = imresize(I,[255 255]) ;
            [im_h, im_w] = size(I);
%               I = double(I);
    
%             if max(im_h, im_w) > maxImSize,
%                 I = imresize(I, maxImSize/max(im_h, im_w), 'bicubic');
%                 [im_h, im_w] = size(I);
%             end;
%             [part part_fea] = fenkuai(I,50,50,50,50);
% histogram=[];
% 
% for i=1:size(part,2) 
%  IPart = part{1,i};
%     
% gaborArray = gaborFilterBank(5,8,39,39);  % Generates the Gabor filter bank
% featureVector = gaborFeatures(I,gaborArray,5,5);



% %     for j=1:3
% %          xx=histeq(IPart(:,:,j));
% %         x=WLDV( double(rgb2gray(IPart)));
% %              siftArr=lbp(rgb2gray(IPart));
% %           y=hist(double(x(:)));
%      siftArr=LBPV(I);
%          histogram=[histogram x];
% %     end
% end

%    siftArr=lpq(I);
%    texturefilters=;
%    filename=['./texturefilters/ICAtextureFilters_11x11_8bit'];
% load(filename, 'ICAtextureFilters');
% texturefilters='E:\Matlab19\ScSPM\texturefilters\ICAtextureFilters_11x11_8bit';
%     mapping=getmapping(16,'u2');
% bins = 2^neighbors;
%    if isstruct(mapping)
% %     bins = mapping.num;
%     for i = 1:size(result,1)
%         for j = 1:size(result,2)
%             result(i,j) = mapping.table(result(i,j)+1);
%         end
%     end
%   end


% 
%  fc = [100;100;10;-pi/8] ;
% % %  [f,d] = vl_sift(I,'Frames', fc);
% % % perm = randperm(size(fc,2)) ;
% % % sel = perm(1:50) ;
% % % h1 = vl_plotframe(fc(:,sel)) ;
% % % h2 = vl_plotframe(fc(:,sel)) ;
%     [f,d] = vl_sift(I,'frames',fc) ;
%        I_       = vl_imsmooth(im2double(I), sqrt(f(3)^2 - 0.5^2)) ;
%           [Ix, Iy] = vl_grad(I_) ;
%           mod      = sqrt(Ix.^2 + Iy.^2) ;
%           ang      = atan2(Iy,Ix) ;
%           grd      = shiftdim(cat(3,mod,ang),2) ;
%           grd      = single(grd) ;
%           siftArr      = vl_siftdescriptor(grd, f) ;
% % 




%  [x,y] =  get_interest_points(I,50);
% *********for PHOG************ % 
%  [y1,x1] = size(I);
% roi = [1;y1;1;x1];
% siftArr=anna_phog(I,20,360,3,roi);
% % siftArr = vl_covdet(patch) ;
% 
% % % % % % % %   siftArr=TPLBP(I);
% siftArr=vl_liop(I);
%  siftArr=LBPTOP(I,[y1,x1],1,3,3,8,3,3,0);
%  siftArr= desc_LMP(I);
% siftArr= calc_LESH(I);
%  corners=detectMSERFeatures(I);
%  [features, valid_corners] = extractFeatures(I, corners);
%  
%   siftArr=features;
figure;
imshow(I)
siftArr=adpmedian(I,3);
figure;
imshow(siftArr)

% % siftArrr=detectKAZEFeatures(I);
% % [siftArr, valid_corners] = extractFeatures(I,siftArrr);
% figure;
% imshow(siftArr);
% 
% % %   siftArr=ldp(I, [y1, x1],10);
%  siftArr=desc_Ltrp(I);
%  figure;
% imshow(siftArr);

%  mapping=getmapping(I,'u2'); 
%       H1=LBP(I,1,8,mapping,'h');
%       siftArr=FPLBP(I);
%  mapping=getmapping(8,'u2'); 
%       H1=LBPV(I,1,8,mapping);
%    siftArr=LBPV(I);
% roi = [1;225;1;300];
% siftArr=anna_phog(I,8,360,3,roi);
% dzy= compute_daisy(I);
% dsc=compute_descriptor(dzy, y, x, ori, spatial_int, hist_int, nt)
%     siftArr=hog_feature_vector(I);
% bins=256;
%   result =mymethod(double(I));
%  siftArr=hist(double(result(:)));
%         siftArr=lpqv(I);
%      siftArr=LBPV(double(I));
%         siftArr=lbp(I);
% siftArr=WLDhist(double(I));

%         siftArr=WLDV(double(I));
%       siftArr=WLDhist(double(I));
%      siftArr=lbp(I);
% siftArr=hist(double(siftArr(:)));
%      siftArr=hist(double(siftArr(:)));
%              feaSet.feaArr = histogram;

 % make grid sampling SIFT descriptors
%             remX = mod(im_w-patchSize,gridSpacing);
%             offsetX = floor(remX/2)+1;
%             remY = mod(im_h-patchSize,gridSpacing);
%             offsetY = floor(remY/2)+1;
%     
%             [gridX,gridY] = meshgrid(offsetX:gridSpacing:im_w-patchSize+1, offsetY:gridSpacing:im_h-patchSize+1);
% 
%             fprintf('Processing %s: wid %d, hgt %d, grid size: %d x %d, %d patches\n', ...
%                      frames(jj).name, im_w, im_h, size(gridX, 2), size(gridX, 1), numel(gridX));
% 
%             % find SIFT descriptors
%             siftArr = sp_find_sift_grid(I, gridX, gridY, patchSize, 0.8); nrml_threshold);
%             
%             siftLens = [siftLens; siftlen];

%                siftArr=WHOG(I);
%      feaSet.feaArr=permute(siftArr,[2 1 3])
          feaSet.feaArr = siftArr';
            feaSet.width = im_w;
            feaSet.height = im_h;
            
%             [siftArr, siftlen] = sp_normalize_sift(siftArr,
            [pdir, fname] = fileparts(frames(jj).name);                        
            fpath = fullfile(rt_data_dir, subname, [fname, '.mat']);
            
            save(fpath, 'feaSet');
            database.path = [database.path, fpath];
            
% 
% % inputFolder = 'E:\Matlab19\ScSPM\data\aa\1';
% outputFolder = 'E:\Matlab19\ScSPM\data'; % Please change, if needed.
%   [pdir, fname] = fileparts(frames(jj).name);   
% fileList = dir(fullfile(rt_data_dir, subname, [fname, '.mat']));
% for kk = 1:numel(fileList)
%   S = load(fullfile(rt_data_dir, subname, fileList(kk).name));
%   I = S.feaSet.feaArr;
%   I = mat2gray(I);
%   fileName = replace(fileList(kk).name,'.mat','.tiff');
%   imwrite(I,fullfile(outputFolder,fileName));
% end

            
%             fid=fopen(fullfile(rt_data_dir, subname, [fname, '.csv']));
%             values= fread(fid);
%   csvwrite(fullfile(rt_data_dir, subname, [fname, '.csv']), values);
%              fclose(fid);
        end;    
    end;
end;
    
lenStat = hist(siftLens, 100);


    
lenStat = hist(siftLens, 100);
