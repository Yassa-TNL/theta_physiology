function [OFC_chan_idx,fro_chan_idx,temp_chan_idx,insula_chan_idx,cingulate_chan_idx,...
    ACC_chan_idx,EC_chan_idx, HC_chan_idx ,CA3_chan_idx,CA1_chan_idx, MTL_chan_idx, NC_chan_idx]  = get_elecs_clean_cue_resp(subj)

cd(['/mnt//yassamri/iEEG/sandra/subj_' subj])
load('sig_elecs_per_reg.mat')
% load cue resp
if strcmp('39', subj)
    OFC_chan_idx       = [1:2, 58:59];
    fro_chan_idx       = [16:19, 25:26, 28, 64, 66, 70:74];
    temp_chan_idx      = [36:37, 48, 57, 93, 102:103];
    insula_chan_idx    = [20];
    cingulate_chan_idx = [11:12];
    ACC_chan_idx       = [11:12];
    CA3_chan_idx       = [40, 87, 95];
    CA1_chan_idx       = [39, 41:42, 52, 86, 94, 96];
    HC_chan_idx        = [39:42, 51, 52, 86:87, 94:96];
    EC_chan_idx        = [76, 85];
    MTL_chan_idx       = [29:32, 39:42, 51, 52, 77:79, 86:87, 94:96];
   
    [length(OFC_chan_idx) length(fro_chan_idx) length(temp_chan_idx)  length(insula_chan_idx)...
  length(cingulate_chan_idx)  length(EC_chan_idx) length(CA3_chan_idx)  length(HC_chan_idx)];

elseif strcmp('44', subj)
    OFC_chan_idx       = [39:40, 77:78, 80];
    fro_chan_idx       = [36, 37, 84:86];
    temp_chan_idx      = [55, 64 67, 74]; % EXCLUDE - BAD ARTIFACT 10, 16, 19, 26:28, 65 67
    insula_chan_idx    = [];
    cingulate_chan_idx = [29:32];
    ACC_chan_idx       = [29:32];
    EC_chan_idx        = [1, 49, 58, 59, 68];
    CA3_chan_idx       = [];
    CA1_chan_idx       = [2, 14, 23];
    HC_chan_idx        = [2, 13, 14, 22, 23];
    MTL_chan_idx       = [2:4, 13:14, 20, 22, 23, 50:51];
   [length(OFC_chan_idx) length(fro_chan_idx) length(temp_chan_idx)  length(insula_chan_idx)...
  length(cingulate_chan_idx)  length(EC_chan_idx) length(CA3_chan_idx)  length(HC_chan_idx)];

elseif strcmp('57', subj)
    OFC_chan_idx       = [18, 19, 85, 87];
    fro_chan_idx       = [1, 7:8, 16:17, 27, 35, 37, 41, 44, 46:47, 79, 81, 83:84, 103, 113 ];
    temp_chan_idx      = [56, 65:66, 73:74, 120, 122, 132, 138:141];
    insula_chan_idx    = [29:32, 38:39, 94:97, 104:105, 118 ];
    cingulate_chan_idx = [2, 9:10, 76:78];
    ACC_chan_idx       = [2];
    EC_chan_idx        = [123];
    CA3_chan_idx       = [59, 60, 133, 134];
    CA1_chan_idx       = [61, 67:68, 125];
    HC_chan_idx        = [59:61, 67:68, 124:125, 133:134];
    MTL_chan_idx       = [48:50, 59:61, 67:68, 114:116, 124:125, 133:134];
   [length(OFC_chan_idx) length(fro_chan_idx) length(temp_chan_idx)  length(insula_chan_idx)...
  length(cingulate_chan_idx)  length(EC_chan_idx) length(CA3_chan_idx)  length(HC_chan_idx)];

elseif strcmp('63', subj)
    OFC_chan_idx       = [37:40, 108:110];
    fro_chan_idx       = [45, 52:56, 116];
    temp_chan_idx      = [10, 14, 23:24, 85:86, 96, 117:120];
    insula_chan_idx    = [];
    cingulate_chan_idx = [47, 98:102];
    ACC_chan_idx       = [47, 98:102];
    EC_chan_idx        = [5:7, 15:16, 67, 91];
    CA3_chan_idx       = [26, 69, 78];
    CA1_chan_idx       = [18, 27, 70, 77, 80, 94];
    HC_chan_idx        = [17:18, 26:27, 68:70, 77:80, 93:94];
    MTL_chan_idx       = [8, 17:18, 25:27, 68:70, 77:80, 93:94];
      [length(OFC_chan_idx) length(fro_chan_idx) length(temp_chan_idx)  length(insula_chan_idx)...
  length(cingulate_chan_idx)  length(EC_chan_idx) length(CA3_chan_idx)  length(HC_chan_idx)];

elseif strcmp('66', subj)
    OFC_chan_idx       = [41, 101];
    fro_chan_idx       = [34, 38:40, 48:50, 58:60, 95:97, 108:110, 118:120];
    temp_chan_idx      = [8:10, 18:20, 68, 70, 78:80];
    insula_chan_idx    = [];
    cingulate_chan_idx = [31:33, 51:52, 91, 111:113];
    ACC_chan_idx       = [31:33, 91];
    EC_chan_idx        = [];
    CA3_chan_idx       = [13, 71:72];
    CA1_chan_idx       = [14, 21:23, 74, 83 ];
    HC_chan_idx        = [12:14, 21:23, 71:74, 81:83];
    MTL_chan_idx       = [1:4, 12:14, 21:23, 61:64, 71:74, 81:83, 85];
      [length(OFC_chan_idx) length(fro_chan_idx) length(temp_chan_idx)  length(insula_chan_idx)...
  length(cingulate_chan_idx)  length(EC_chan_idx) length(CA3_chan_idx)  length(HC_chan_idx)];

elseif strcmp('83', subj)
    OFC_chan_idx       = [];
    fro_chan_idx       = [8:12, 15:20, 23:28, 32:36, 40:43, 49:50, ...
        57:58, 75:79, 91:95];
    temp_chan_idx      = [5:6, 13, 21];
    insula_chan_idx    = [];
    cingulate_chan_idx = [67:68, 71:74, 83:84, 87:90];
    ACC_chan_idx       = [67:68, 83:84];
    EC_chan_idx        = [];
    CA3_chan_idx       = [];
    CA1_chan_idx       = [];
    HC_chan_idx        = [];
    MTL_chan_idx       = [];
      [length(OFC_chan_idx) length(fro_chan_idx) length(temp_chan_idx)  length(insula_chan_idx)...
  length(cingulate_chan_idx)  length(EC_chan_idx) length(CA3_chan_idx)  length(HC_chan_idx)];

elseif strcmp('84', subj) % this section is all wrong. Dont know whats goin on here.
    OFC_chan_idx       = [10, 48, 50, 51];
    fro_chan_idx       = [36, 47, 74, 81, 87:90, 107, 110:111];
    temp_chan_idx      = [33, 68, 70, 98:101, 140:141, 143, 153];
    insula_chan_idx    = [104, 114:115];
    cingulate_chan_idx = [82];
    ACC_chan_idx       = [];
    EC_chan_idx        = [64:65, 146];
    CA3_chan_idx       = [26];
    CA1_chan_idx       = [27:28, 92, 134, 144];
    HC_chan_idx        = [26:28, 92, 134, 144];
    MTL_chan_idx       = [15:17, 26:28, 92, 56:59, 66, 134, 144];
        [length(OFC_chan_idx) length(fro_chan_idx) length(temp_chan_idx)  length(insula_chan_idx)...
  length(cingulate_chan_idx)  length(EC_chan_idx) length(CA3_chan_idx)  length(HC_chan_idx)] 

elseif strcmp('85', subj)
    OFC_chan_idx       = [48:49];
    fro_chan_idx       = [55, 70, 91, 93];
    temp_chan_idx      = [5, 7:8, 14, 22, 33, 82:83, 96:97, 99:103, 110:111];
    insula_chan_idx    = [66:67];
    cingulate_chan_idx = [34, 56, 57, 84];
    ACC_chan_idx       = [34];
    EC_chan_idx        = [77, 106];
    CA3_chan_idx       = [];
    CA1_chan_idx       = [28, 105];
    HC_chan_idx        = [26:28, 75, 104, 105];
    MTL_chan_idx       = [1:2, 9:10, 15:17, 23, 26:28, 75, 104, 105];
          [length(OFC_chan_idx) length(fro_chan_idx) length(temp_chan_idx)  length(insula_chan_idx)...
  length(cingulate_chan_idx)  length(EC_chan_idx) length(CA3_chan_idx)  length(HC_chan_idx)];   

elseif strcmp('87', subj)
    OFC_chan_idx       = [64];
    fro_chan_idx       = [22, 25, 59, 61:63, 88:91];
    temp_chan_idx      = [12:13, 43, 46:47, 52:54, 79, 81, 96, 97, 100];
    insula_chan_idx    = [82, 84, 85, 103:105];
    cingulate_chan_idx = [16:17, 57];
    ACC_chan_idx       = [16:17, 57];
    EC_chan_idx        = [];
    CA3_chan_idx       = [];
    CA1_chan_idx       = [93];
    HC_chan_idx        = [9, 72, 73, 93];
    MTL_chan_idx       = [1:2, 9, 34:35, 72:74, 93];
              [length(OFC_chan_idx) length(fro_chan_idx) length(temp_chan_idx)  length(insula_chan_idx)...
  length(cingulate_chan_idx)  length(EC_chan_idx) length(CA3_chan_idx)  length(HC_chan_idx)];   

end
    OFC_chan_idx       = OFC_chan_idx(logical(sig_elecs_per_reg{1}));
    fro_chan_idx       = fro_chan_idx(logical(sig_elecs_per_reg{2}));
    temp_chan_idx      = temp_chan_idx(logical(sig_elecs_per_reg{3}));
    insula_chan_idx    = insula_chan_idx(logical(sig_elecs_per_reg{4}));
    cingulate_chan_idx = cingulate_chan_idx(logical(sig_elecs_per_reg{5}));
    EC_chan_idx        = EC_chan_idx(logical(sig_elecs_per_reg{6}));
    HC_chan_idx        = HC_chan_idx(logical(sig_elecs_per_reg{7}));
    CA3_chan_idx       = CA3_chan_idx(logical(sig_elecs_per_reg{8}));

    NC_chan_idx        = [OFC_chan_idx fro_chan_idx temp_chan_idx insula_chan_idx cingulate_chan_idx EC_chan_idx]; 

        [length(OFC_chan_idx) length(fro_chan_idx) length(temp_chan_idx)  length(insula_chan_idx)...
  length(cingulate_chan_idx)  length(EC_chan_idx) length(CA3_chan_idx)  length(HC_chan_idx)];
%%

% if strcmp('39', subj)
%     OFC_chan_idx       = [1:2];
%     fro_chan_idx       = [64    72    74];
%     temp_chan_idx      = [49 93 103];
%     insula_chan_idx    = [];
%     cingulate_chan_idx = [];
%     ACC_chan_idx       = [];
%     CA3_chan_idx       = [40, 87, 95];
%     CA1_chan_idx       = [];
%     HC_chan_idx        = [40    42    86    87];
%     EC_chan_idx        = [];
%     MTL_chan_idx       = [];
%     [length(OFC_chan_idx) length(fro_chan_idx) length(temp_chan_idx)  length(insula_chan_idx)...
%     length(cingulate_chan_idx)  length(EC_chan_idx) length(CA3_chan_idx)  length(HC_chan_idx)];
% 
% elseif strcmp('44', subj)
%     OFC_chan_idx       = [ ];
%     fro_chan_idx       = [85:86];
%     temp_chan_idx      = []; 
%     insula_chan_idx    = [];
%     cingulate_chan_idx = [31];
%     ACC_chan_idx       = [];
%     EC_chan_idx        = [58, 68];
%     CA3_chan_idx       = [];
%     CA1_chan_idx       = [];
%     HC_chan_idx        = [];
%     MTL_chan_idx       = [];
%    [length(OFC_chan_idx) length(fro_chan_idx) length(temp_chan_idx)  length(insula_chan_idx)...
%    length(cingulate_chan_idx)  length(EC_chan_idx) length(CA3_chan_idx)  length(HC_chan_idx)];
% 
% elseif strcmp('57', subj)
%     OFC_chan_idx       = [18, 19, 87];
%     fro_chan_idx       = [1, 7:8, 16:17,  35,  44, 46:47, 79, 81, 83:84, 103];
%     temp_chan_idx      = [56, 65, 73:74, 122, 132];
%     insula_chan_idx    = [31    38    39    97   104   105   118];
%     cingulate_chan_idx = [2, 9:10];
%     ACC_chan_idx       = [];
%     EC_chan_idx        = [];
%     CA3_chan_idx       = [59];
%     CA1_chan_idx       = [];
%     HC_chan_idx        = [59, 67:68];
%     MTL_chan_idx       = [];
%    [length(OFC_chan_idx) length(fro_chan_idx) length(temp_chan_idx)  length(insula_chan_idx)...
%    length(cingulate_chan_idx)  length(EC_chan_idx) length(CA3_chan_idx)  length(HC_chan_idx)];
% 
% elseif strcmp('63', subj)
%     OFC_chan_idx       = [37:40, 110];
%     fro_chan_idx       = [45, 52:56, 116];
%     temp_chan_idx      = [10,  23:24, 86, 96, 117:120];
%     insula_chan_idx    = [];
%     cingulate_chan_idx = [99];
%     ACC_chan_idx       = [];
%     EC_chan_idx        = [5:6, 15:16];
%     CA3_chan_idx       = [78];
%     CA1_chan_idx       = [];
%     HC_chan_idx        = [17    77    78    94];
%     MTL_chan_idx       = [];
%       [length(OFC_chan_idx) length(fro_chan_idx) length(temp_chan_idx)  length(insula_chan_idx)...
%   length(cingulate_chan_idx)  length(EC_chan_idx) length(CA3_chan_idx)  length(HC_chan_idx)];
% 
% elseif strcmp('66', subj)
%     OFC_chan_idx       = [41, 101];
%     fro_chan_idx       = [34, 38:40, 48:50, 59    96   108   109   110 118:120];
%     temp_chan_idx      = [8:10, 18:20, 68, 70, 80];
%     insula_chan_idx    = [];
%     cingulate_chan_idx = [31:33, 51:52, 91, 111];
%     ACC_chan_idx       = [];
%     EC_chan_idx        = [];
%     CA3_chan_idx       = [13];
%     CA1_chan_idx       = [];
%     HC_chan_idx        = [13:14, 21:22, 73, 82];
%     MTL_chan_idx       = [];
%       [length(OFC_chan_idx) length(fro_chan_idx) length(temp_chan_idx)  length(insula_chan_idx)...
%   length(cingulate_chan_idx)  length(EC_chan_idx) length(CA3_chan_idx)  length(HC_chan_idx)];
% 
% elseif strcmp('83', subj)
%     OFC_chan_idx       = [];
%     fro_chan_idx       = [];
%     temp_chan_idx      = [];
%     insula_chan_idx    = [];
%     cingulate_chan_idx = [];
%     ACC_chan_idx       = [];
%     EC_chan_idx        = [];
%     CA3_chan_idx       = [];
%     CA1_chan_idx       = [];
%     HC_chan_idx        = [];
%     MTL_chan_idx       = [];
%       [length(OFC_chan_idx) length(fro_chan_idx) length(temp_chan_idx)  length(insula_chan_idx)...
%   length(cingulate_chan_idx)  length(EC_chan_idx) length(CA3_chan_idx)  length(HC_chan_idx)];
% 
% 
% elseif strcmp('84', subj) 
%     OFC_chan_idx       = [10, 48];
%     fro_chan_idx       = [36, 87, 107];
%     temp_chan_idx      = [70, 143, 153];
%     insula_chan_idx    = [115];
%     cingulate_chan_idx = [82];
%     ACC_chan_idx       = [];
%     EC_chan_idx        = [146];
%     CA3_chan_idx       = [];
%     CA1_chan_idx       = [];
%     HC_chan_idx        = [27:28, 134, 144];
%     MTL_chan_idx       = [];
%         [length(OFC_chan_idx) length(fro_chan_idx) length(temp_chan_idx)  length(insula_chan_idx)...
%   length(cingulate_chan_idx)  length(EC_chan_idx) length(CA3_chan_idx)  length(HC_chan_idx)] 
% 
% 
% 
% elseif strcmp('85', subj)
%     OFC_chan_idx       = [48:49];
%     fro_chan_idx       = [55, 70, 91];
%     temp_chan_idx      = [33    82    83    96];
%     insula_chan_idx    = [66];
%     cingulate_chan_idx = [57];
%     ACC_chan_idx       = [34];
%     EC_chan_idx        = [106];
%     CA3_chan_idx       = [];
%     CA1_chan_idx       = [];
%     HC_chan_idx        = [26 75, 104, 105];
%     MTL_chan_idx       = [];
%           [length(OFC_chan_idx) length(fro_chan_idx) length(temp_chan_idx)  length(insula_chan_idx)...
%   length(cingulate_chan_idx)  length(EC_chan_idx) length(CA3_chan_idx)  length(HC_chan_idx)]; 
% 
% 
% elseif strcmp('87', subj)
%     OFC_chan_idx       = [64];
%     fro_chan_idx       = [22, 25, 59, 61:63, 88:91];
%     temp_chan_idx      = [12:13, 46:47, 52:54, 79, 81];
%     insula_chan_idx    = [82, 85, 105];
%     cingulate_chan_idx = [17, 57];
%     ACC_chan_idx       = [];
%     EC_chan_idx        = [];
%     CA3_chan_idx       = [];
%     CA1_chan_idx       = [93];
%     HC_chan_idx        = [9, 72, 73, 93];
%     MTL_chan_idx       = [];
%               [length(OFC_chan_idx) length(fro_chan_idx) length(temp_chan_idx)  length(insula_chan_idx)...
%   length(cingulate_chan_idx)  length(EC_chan_idx) length(CA3_chan_idx)  length(HC_chan_idx)];   


% 
% end
% 
%  NC_chan_idx        = [OFC_chan_idx fro_chan_idx temp_chan_idx insula_chan_idx cingulate_chan_idx EC_chan_idx]; 
% 

end

