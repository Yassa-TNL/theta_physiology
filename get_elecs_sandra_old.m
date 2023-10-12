function [OFC_chan_idx,fro_chan_idx,temp_chan_idx,insula_chan_idx,cingulate_chan_idx,...
    ACC_chan_idx,EC_chan_idx, HC_chan_idx ,CA3_chan_idx,CA1_chan_idx, MTL_chan_idx]  = get_elecs_clean(subj)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if strcmp('57', subj)
    fro_chan_idx       = [1 4:8 13 16:17 25:27 35 37 44 46:47 70:71 78:83 89 91 98 99 101 102 108 109 110 112];
    MTL_chan_idx       = [48:50 58:61 67:68 113:115 122:124 132:134];
    temp_chan_idx      = [55 56 65 66 72 73 118:121 128:131 136:140];
    insula_chan_idx    = [29:32 38:40 93:96 103 104];
    cingulate_chan_idx = [2 9:11 74:77];
    OFC_chan_idx       = [18:19 84:86];
    CA3_chan_idx       = [59 60 62 124 132 133];
    CA1_chan_idx       = [58 67 68  134];
    HC_chan_idx        = [59 60 62 124 132 133 134 58 67 68];
    ACC_chan_idx       = [2];
    EC_chan_idx = [];
elseif strcmp('39', subj)
    fro_chan_idx       = [8 16:19 25:26 28 63:64 66 72:74];
    MTL_chan_idx       = [29:32 39:42 51:52 77:79 85:88 94:96];
    temp_chan_idx      = [34:38 45:48 54:57 83:84 91:93 100:103 ];
    insula_chan_idx    = [20 21];
    cingulate_chan_idx = [11 12];
    OFC_chan_idx       = [1 2 58];
    CA3_chan_idx       = [40 41 51 87 88 95];
    CA1_chan_idx       = [39 42 52 85 86 94 96];
    HC_chan_idx        = [40 41 51 87 88 95 39 42 52 85 86 94 96];
    ACC_chan_idx       = [11 12];
    EC_chan_idx = [];
elseif strcmp('66', subj)
    fro_chan_idx       = [37:40 49:50 95:97 101 108 109 110];
    MTL_chan_idx       = [1:5 11:15 21:23 61:65 71:75 81:83 ];
    temp_chan_idx      = [7:9 16:20 24 25 29 30 48 68:70 77 78 84:90 ];
    insula_chan_idx    = [];
    cingulate_chan_idx = [31:33 35 51:54 91:93 111:114];
    OFC_chan_idx       = [41 42 102];
    CA3_chan_idx       = [ 13 23 82 ];
    CA1_chan_idx       = [12 14 22 81 83];
    HC_chan_idx        = [13 23 82 12 14 22 81 83];
    ACC_chan_idx       = [31 32 33 35];
    EC_chan_idx = [];
elseif strcmp('63', subj)
    fro_chan_idx       = [45 52:56 116 ];
    MTL_chan_idx       = [5:8 17:19 25:28 67:70 77:80 91:94];
    temp_chan_idx      = [14 22 23 24 31 34 73:76 84:86];
    insula_chan_idx    = [];
    cingulate_chan_idx = [47 49 97:101];
    OFC_chan_idx       = [37:39];
    CA3_chan_idx       = [  ];
    CA1_chan_idx       = [  ];
    HC_chan_idx        = [];
    ACC_chan_idx      = [47 49];
    EC_chan_idx = [];
elseif strcmp('44', subj)
    fro_chan_idx       = [36:38 44:47 84 :86];
    MTL_chan_idx       = [1:4 11:14 21 22 23 49 58:62 68 70:72];
    temp_chan_idx      = [7:10 16:19 25:28 54:57 64:67 74:76];
    insula_chan_idx    = [];
    cingulate_chan_idx = [29:32];
    OFC_chan_idx       = [39 40 77:80];
    CA3_chan_idx       = [22  ];
    CA1_chan_idx       = [2 3 12 13 14 23 61 62 70 72];
    HC_chan_idx        = [22 2 3 12 13 14 23 61 62 70 72];
    ACC_chan_idx       = [29:32];
    EC_chan_idx = [];
elseif strcmp('83', subj)
    fro_chan_idx       = [7:12 15:20 24:28 33:36 41:44 49:52 59:60 65:66];
    MTL_chan_idx       = [];
    temp_chan_idx      = [5 6 13 21];
    insula_chan_idx    = [];
    cingulate_chan_idx = [];
    OFC_chan_idx       = [];
    CA3_chan_idx = [];
    CA1_chan_idx = [];
    HC_chan_idx        = [];
    ACC_chan_idx       = [];
elseif strcmp('84', subj)
    fro_chan_idx       = [27 29 31 32  28 30 ];
    MTL_chan_idx       = [17:20 41:44 49:52];
    temp_chan_idx      = [11:13 15:16 44 47 48 55 56];
    insula_chan_idx    = [];
    cingulate_chan_idx = [];
    OFC_chan_idx       = [3 4 33 34 35 40];
    CA3_chan_idx       = [ ];
    CA1_chan_idx       = [41:43];
    HC_chan_idx        = [41:43];
    ACC_chan_idx       = [];
elseif strcmp('85', subj)
    fro_chan_idx       = [23:24 46:48];
    MTL_chan_idx       = [1:2 9:11 25:27 33:35];
    temp_chan_idx      = [5:8 14:16 31:32 40];
    insula_chan_idx    = [4];
    cingulate_chan_idx = [17:18];
    OFC_chan_idx       = [41:43];
    CA3_chan_idx       = [];
    CA1_chan_idx       = [35];
    HC_chan_idx        = [35];
    ACC_chan_idx       = [17:18];
elseif strcmp('87', subj)
    fro_chan_idx       = [];
    MTL_chan_idx       = [];
    temp_chan_idx      = [];
    insula_chan_idx    = [];
    cingulate_chan_idx = [];
    OFC_chan_idx       = [];
    CA3_chan_idx       = [];
    CA1_chan_idx       = [];
    HC_chan_idx        = [];
    ACC_chan_idx       = [];
end
end

