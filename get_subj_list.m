function [subj_list] = get_subj_list(reg)
% output: list of subject w/ elecs in reg
if strcmp('OFC',reg) || strcmp('MTL',reg)
    subj_list = {'39' '57' '44' '63' '66' '85' '87' '84'}
elseif strcmp('FRO',reg) || strcmp('TEMP',reg) || strcmp('CING',reg)
    subj_list = {'39' '57' '44' '63' '66'  '85' '87' '84'}
elseif strcmp('ins',reg)
    subj_list = {'39' '57' '84' '85' '87'}
elseif strcmp('EC',reg)
    subj_list = {'39' '44' '57' '63' '84' '85'}
elseif strcmp('CA3',reg)
    subj_list = {'39' '57'  '63' '66' '84' }
elseif strcmp('CA1',reg) || strcmp('HC',reg) || strcmp('NC',reg)
    subj_list = {'39' '44' '57' '63' '66' '84' '85' '87'}
elseif strcmp('OFC_FRO_TEMP',reg)
    subj_list = {'39' '44' '57' '63' '66' '84' '85' '87'}    
end
end

