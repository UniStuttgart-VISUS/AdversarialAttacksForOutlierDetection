function corr_coef=fisher_z_transformation(V)
          for i=1:size(V,2)
             current_corr=V(setdiff(1:size(V,2),i),i);
             z=nanmean(current_corr);
             corr_coef(i)=tanh(z);     
          end
end