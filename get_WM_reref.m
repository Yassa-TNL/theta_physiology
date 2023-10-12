function [WM_ref_vector] = get_WM_reref(subj, chan_label)

WM_ref_vector = zeros(length(chan_label),1);

if strcmp('39', subj)
    WM_ref_vector(1:10)  = 5;
    WM_ref_vector(11:19) = 14;
    WM_ref_vector(20:28) = 23;
    WM_ref_vector(29:38) = 33;
    WM_ref_vector(39:48) = 43;
    WM_ref_vector(49:57) = 53;
    WM_ref_vector(58:66) = 62;
    WM_ref_vector(67:75) = 68;
    WM_ref_vector(76:84) = 80;
    WM_ref_vector(85:93) = 89;
    WM_ref_vector(94:103)= 98;
    
elseif strcmp('44', subj)
    WM_ref_vector(1:10)  = 6;
    WM_ref_vector(11:19) = 15;
    WM_ref_vector(20:28) = 24;
    WM_ref_vector(29:38) = 33;
    WM_ref_vector(39:48) = 41;
    WM_ref_vector(49:57) = 53;
    WM_ref_vector(58:67) = 63;
    WM_ref_vector(68:76) = 73;
    WM_ref_vector(77:86) = 81;
    
elseif strcmp('57', subj)
    WM_ref_vector(1:8)     =3;
    WM_ref_vector(9:17)    =12;
    WM_ref_vector(18:27)   =23;
    WM_ref_vector(28:37)   =33;
    WM_ref_vector(38:47)   =42;
    WM_ref_vector(48:57)   =52;
    WM_ref_vector(58:66)   =63;
    WM_ref_vector(67:74)   =70;
    WM_ref_vector(75:84)   =80;
    WM_ref_vector(85:93)   =88;
    WM_ref_vector(94:103)  =98;
    WM_ref_vector(104:113) =107;
    WM_ref_vector(114:122) =107;
    WM_ref_vector(123:132) =128;
    WM_ref_vector(133:141) =136;
    
elseif strcmp('63', subj)
    WM_ref_vector(1:14)    =12;
    WM_ref_vector(15:24)   =21;
    WM_ref_vector(25:36)   =29;
    WM_ref_vector(37:46)   =42;
    WM_ref_vector(47:66)   =51;
    WM_ref_vector(67:76)   =73;
    WM_ref_vector(77:86)   =82;
    WM_ref_vector(87:96)   =82;
    WM_ref_vector(97:106)  =104;
    WM_ref_vector(107:116) =112;
    WM_ref_vector(117:127) =124;
    
elseif strcmp('66', subj)
    WM_ref_vector(1:10)    =6;
    WM_ref_vector(11:20)   =6;
    WM_ref_vector(21:30)   =27;
    WM_ref_vector(31:40)   =36;
    WM_ref_vector(41:50)   =45;
    WM_ref_vector(51:62)   =55;
    WM_ref_vector(63:70)   =66;
    WM_ref_vector(71:80)   =87;
    WM_ref_vector(81:90)   =87;
    WM_ref_vector(91:100)  =93;
    WM_ref_vector(101:110) =106;
    WM_ref_vector(111:120) =116;
    
    % elseif strcmp('83', subj)
    
elseif strcmp('84', subj)
    WM_ref_vector(1:8) =2;
    WM_ref_vector(9:14) =11;
    WM_ref_vector(15:22) =19;
    WM_ref_vector(23:24) =11;
    WM_ref_vector(25) = 19;
    WM_ref_vector(26:33) =30;
    WM_ref_vector(34:47) =34;
    WM_ref_vector(48:55) =54;
    WM_ref_vector(56:63) =61;
    WM_ref_vector(64:71) =69;
    WM_ref_vector(72:81) =77;
    WM_ref_vector(82:91) =85;
    WM_ref_vector(92:101) =94;
    WM_ref_vector(102:111) =106;
    WM_ref_vector(112:133) =112;
    WM_ref_vector(134:143) =138;
    WM_ref_vector(144:154) =147;

elseif strcmp('85', subj)
    WM_ref_vector(1:8) =3;
    WM_ref_vector(9:14) =11;
    WM_ref_vector(15:22) =19;
    WM_ref_vector(23:33) =31;
    WM_ref_vector(34:47) =36;
    WM_ref_vector(48:55) =51;
    WM_ref_vector(56:63) =60;
    WM_ref_vector(64:73) =68;
    WM_ref_vector(74:83) =80;
    WM_ref_vector(84:93) =88;
    WM_ref_vector(94:103) =80;
    WM_ref_vector(104:126) =108;
    
elseif strcmp('87', subj)
    WM_ref_vector(1:8) =5;
    WM_ref_vector(9:14) =11;
    WM_ref_vector(15:22) =18;
    WM_ref_vector(23:24) =11;
    WM_ref_vector(25) =18;
    WM_ref_vector(26:33) =30;
    WM_ref_vector(34:47) =37;
    WM_ref_vector(48:55) =50;
    WM_ref_vector(56:63) =60;
    WM_ref_vector(64:71) =69;
    WM_ref_vector(72:81) =77;
    WM_ref_vector(82:91) =86;
    WM_ref_vector(92:101) =98;
    WM_ref_vector(102:124) =107;
    
end
end


