%% my own code written up for pte paired non parametric perm testing   

if iCond ==1
    real_diff = mean(reg1_reg2_condc - reg1_reg2_condb);
    both_conds = cat(1, reg1_reg2_condc, reg1_reg2_condb);
    elseif iCond ==2 
    real_diff = mean(reg1_reg2_conda - reg1_reg2_condb);
    both_conds = cat(1, reg1_reg2_conda, reg1_reg2_condb);        
    elseif iCond ==3
    real_diff = mean(reg1_reg2_condd - reg1_reg2_condb);
    both_conds = cat(1, reg1_reg2_condd, reg1_reg2_condb);        
    end
    null_dist_of_diff = nan(1,n_permutes);
    for permi = 1:n_permutes
        shuff_cntr = mod(reshape(randperm(1*size(both_conds,1)), 1, size(both_conds,1)), 2 );
        cond1 =both_conds(logical(shuff_cntr));
        cond2 =both_conds(logical(shuff_cntr==0));
        null_dist_of_diff (permi) = mean(cond1-cond2);
    end
    
    p1a = length(find(null_dist_of_diff>real_diff))/n_permutes;
    p1b = length(find(null_dist_of_diff<real_diff))/n_permutes;
    
    if real_diff<0
       p_val (iCond) = p1b;
   else
       p_val (iCond) = p1b;
    end