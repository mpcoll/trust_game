X = (0:3);%[0]% 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23];
%Y = (23:47);
Y = 4;
IDisapp = zeros(1,4);
%n = 216875;
%n = 224640;
n=20;
tgrid = [4 3 2 1];
lgrid = [0 2 4 6 8];
k = 93;
IHold = zeros(1,k);
THold = zeros(1,k);
IoHold = zeros(1,k);
ToHold = zeros(1,k);
IavHold = zeros(1,k);
TavHold = zeros(1,k);
ItempHold = zeros(1,k);
TtempHold = zeros(1,k);
IirrHold = zeros(1,k);
TirrHold = zeros(1,k);
subj = 0;
for j = 1:(length(X))
fid =fopen(['VarEstIrr' int2str(X(j)) '.bin'],'r'); %SampleGenerative7
for i = 1:n
    %if i>1
    %    subj = subj+1;
    %end
    f = 0;
    s = 0;
    u =0;
    o=0;
    su = 0;
    w =0;
    g=0;   
    he = 0;
    ue = 0;
    subj = subj+1;
   %Keep =  fread(fid,1,'double');
   for risk = 0:7
      % s = s+1;
      for t = 1:4
         % s = s+1;
          for plan = 1:4
              %s = s+1;
              for shift = 0:4
                 % s = s+1;
                 for level = 1:3
                     %s = s+1;
                     for irr = 0:4
                         %s = s+1;
                         for guilt = 0:2
                             s = s+1;
                             Keep =  fread(fid,1,'double');
                            % totalprob(subj) = totalprob(subj)+exp(-Keep);
                            if exp(-Keep) > IHold(subj)
                                IHold(subj) = exp(-Keep);
                                MInvestorLikelihood(subj)= Keep;
                                MInvestorPlan(subj)=plan;
                                MInvestorToM(subj)=2*(level-1);                             
                                MInvestorAversion(subj) = risk;
                                MInvestorGuilt(subj) = guilt;
                                MInvestorTemp(subj) = tgrid(t);
                                MInvestorIrritability(subj) = irr;
                                MInvestorShift(subj) = shift;
                            end
                             if ( shift == 0)
                                 ue = ue+1;
                                if exp(-Keep) > IirrHold(subj)                                 
                                IirrHold(subj) = exp(-Keep); 
                                IInvestorLikelihood(subj)= Keep;
                                IInvestorPlan(subj)=plan;
                                IInvestorToM(subj)=2*(level-1);                             
                                IInvestorAversion(subj) = risk;
                                IInvestorGuilt(subj) = guilt;
                                IInvestorTemp(subj) = tgrid(t);
                                IInvestorIrritability(subj) = irr;
                                IInvestorShift(subj) = shift;                                
                                end
                             end                               
                             if (irr == 0 & shift == 0)
                                 u = u+1;
                                if exp(-Keep) > IavHold(subj)                                 
                                IavHold(subj) = exp(-Keep); 
                                AvInvestorLikelihood(subj)= Keep;
                                AvInvestorPlan(subj)=plan;
                                AvInvestorToM(subj)=2*(level-1);                             
                                AvInvestorAversion(subj) = risk;
                                AvInvestorGuilt(subj) = guilt;
                                AvInvestorTemp(subj) = tgrid(t);
                                AvInvestorIrritability(subj) = irr;
                                AvInvestorShift(subj) = shift;                                
                                end
                             end      
                             if (irr == 0 & risk == 3 & t== 2 & shift == 0)
                                 o = o+1;
                                if exp(-Keep) > IoHold(subj)                                 
                                IoHold(subj) = exp(-Keep); 
                                OInvestorLikelihood(subj)= Keep;
                                OInvestorPlan(subj)=plan;
                                OInvestorToM(subj)=2*(level-1);                             
                                OInvestorAversion(subj) = risk;
                                OInvestorGuilt(subj) = guilt;
                                OInvestorTemp(subj) = tgrid(t);
                                OInvestorIrritability(subj) = irr;
                                OInvestorShift(subj) = shift;                                
                                end
                             end      
                             if (irr == 0 & risk == 3 & shift == 0)
                                 he = he+1;
                                if exp(-Keep) > ItempHold(subj)                                 
                                ItempHold(subj) = exp(-Keep); 
                                TInvestorLikelihood(subj)= Keep;
                                TInvestorPlan(subj)=plan;
                                TInvestorToM(subj)=2*(level-1);                             
                                TInvestorAversion(subj) = risk;
                                TInvestorGuilt(subj) = guilt;
                                TInvestorTemp(subj) = tgrid(t);
                                TInvestorIrritability(subj) = irr;
                                TInvestorShift(subj) = shift;                                
                                end
                             end                               
                         end
                     end 
                     for irr = 0:4
                   
                        if ne(level,3)
                         for guilt = 0:2
                             Keep =  fread(fid,1,'double');
                             su =  su+1; 
                             if exp(-Keep) > THold(subj)
                                THold(subj) = exp(-Keep);
                                MTrusteeLikelihood(subj)= Keep;
                                MTrusteePlan(subj)=plan;
                                MTrusteeToM(subj)=2*level-1;                             
                                MTrusteeAversion(subj) = risk;
                                MTrusteeGuilt(subj) = guilt;
                                MTrusteeTemp(subj) = tgrid(t);
                                MTrusteeIrritability(subj) = irr;
                                MTrusteeShift(subj) = shift;
                             end  
                              if ( shift == 0)
                                 ue = ue+1;
                                if exp(-Keep) > TirrHold(subj)                                 
                                TirrHold(subj) = exp(-Keep); 
                                ITrusteeLikelihood(subj)= Keep;
                                ITrusteePlan(subj)=plan;
                                ITrusteeToM(subj)=2*(level-1);                             
                                ITrusteeAversion(subj) = risk;
                                ITrusteeGuilt(subj) = guilt;
                                ITrusteeTemp(subj) = tgrid(t);
                                ITrusteeIrritability(subj) = irr;
                                ITrusteeShift(subj) = shift;                                
                                end
                             end                                
                             if (irr == 0 & shift == 0)
                                 u = u+1;
                                if exp(-Keep) > TavHold(subj)                                 
                                TavHold(subj) = exp(-Keep); 
                                AvTrusteeLikelihood(subj)= Keep;
                                AvTrusteePlan(subj)=plan;
                                AvTrusteeToM(subj)=2*(level-1);                             
                                AvTrusteeAversion(subj) = risk;
                                AvTrusteeGuilt(subj) = guilt;
                                AvTrusteeTemp(subj) = tgrid(t);
                                AvTrusteeIrritability(subj) = irr;
                                AvTrusteeShift(subj) = shift;                                
                                end
                             end      
                             if (irr == 0 & risk == 3 & t== 2 & shift == 0)
                                 o = o+1;
                                if exp(-Keep) > ToHold(subj)                                 
                                ToHold(subj) = exp(-Keep); 
                                OTrusteeLikelihood(subj)= Keep;
                                OTrusteePlan(subj)=plan;
                                OTrusteeToM(subj)=2*(level-1);                             
                                OTrusteeAversion(subj) = risk;
                                OTrusteeGuilt(subj) = guilt;
                                OTrusteeTemp(subj) = tgrid(t);
                                OTrusteeIrritability(subj) = irr;
                                OTrusteeShift(subj) = shift;                                
                                end
                             end          
                             if (irr == 0 & risk == 3  & shift == 0)
                                 he = he+1;
                                if exp(-Keep) > TtempHold(subj)                                 
                                TtempHold(subj) = exp(-Keep); 
                                TTrusteeLikelihood(subj)= Keep;
                                TTrusteePlan(subj)=plan;
                                TTrusteeToM(subj)=2*(level-1);                             
                                TTrusteeAversion(subj) = risk;
                                TTrusteeGuilt(subj) = guilt;
                                TTrusteeTemp(subj) = tgrid(t);
                                TTrusteeIrritability(subj) = irr;
                                TTrusteeShift(subj) = shift;                                
                                end
                             end                                  
                         end      
                        end
                     end
                     
                 end


                     
                                
              end
          end
      end
   end
end
    fclose(fid);
end

n =13;
for j = 1:(length(Y))
fid =fopen(['VarEstIrr' int2str(Y(j)) '.bin'],'r'); %SampleGenerative7
for i = 1:n
    %if i>1
    %    subj = subj+1;
    %end
    f = 0;
    s = 0;
    u =0;
    o=0;
    su = 0;
    w =0;
    g=0;    
    he = 0;
    ue = 0;
    subj = subj+1;
   %Keep =  fread(fid,1,'double');
   for risk = 0:7
      % s = s+1;
      for t = 1:4
         % s = s+1;
          for plan = 1:4
              %s = s+1;
              for shift = 0:4
                 % s = s+1;
                 for level = 1:3
                     %s = s+1;
                     for irr = 0:4
                         %s = s+1;
                         for guilt = 0:2
                             s = s+1;
                             Keep =  fread(fid,1,'double');
                            % totalprob(subj) = totalprob(subj)+exp(-Keep);
                            if exp(-Keep) > IHold(subj)
                                IHold(subj) = exp(-Keep);
                                MInvestorLikelihood(subj)= Keep;
                                MInvestorPlan(subj)=plan;
                                MInvestorToM(subj)=2*(level-1);                             
                                MInvestorAversion(subj) = risk;
                                MInvestorGuilt(subj) = guilt;
                                MInvestorTemp(subj) = tgrid(t);
                                MInvestorIrritability(subj) = irr;
                                MInvestorShift(subj) = shift;
                            end
                             if ( shift == 0)
                                 ue = ue+1;
                                if exp(-Keep) > IirrHold(subj)                                 
                                IirrHold(subj) = exp(-Keep); 
                                IInvestorLikelihood(subj)= Keep;
                                IInvestorPlan(subj)=plan;
                                IInvestorToM(subj)=2*(level-1);                             
                                IInvestorAversion(subj) = risk;
                                IInvestorGuilt(subj) = guilt;
                                IInvestorTemp(subj) = tgrid(t);
                                IInvestorIrritability(subj) = irr;
                                IInvestorShift(subj) = shift;                                
                                end
                             end                             
                             if (irr == 0 & shift == 0)
                                 u = u+1;
                                if exp(-Keep) > IavHold(subj)                                 
                                IavHold(subj) = exp(-Keep); 
                                AvInvestorLikelihood(subj)= Keep;
                                AvInvestorPlan(subj)=plan;
                                AvInvestorToM(subj)=2*(level-1);                             
                                AvInvestorAversion(subj) = risk;
                                AvInvestorGuilt(subj) = guilt;
                                AvInvestorTemp(subj) = tgrid(t);
                                AvInvestorIrritability(subj) = irr;
                                AvInvestorShift(subj) = shift;                                
                                end
                             end      
                             if (irr == 0 & risk == 3 & t== 2 & shift == 0)
                                 o = o+1;
                                if exp(-Keep) > IoHold(subj)                                 
                                IoHold(subj) = exp(-Keep); 
                                OInvestorLikelihood(subj)= Keep;
                                OInvestorPlan(subj)=plan;
                                OInvestorToM(subj)=2*(level-1);                             
                                OInvestorAversion(subj) = risk;
                                OInvestorGuilt(subj) = guilt;
                                OInvestorTemp(subj) = tgrid(t);
                                OInvestorIrritability(subj) = irr;
                                OInvestorShift(subj) = shift;                                
                                end
                             end          
                             if (irr == 0 & risk == 3 & t== 2 & shift == 0)
                                 he = he+1;
                                if exp(-Keep) > ItempHold(subj)                                 
                                ItempHold(subj) = exp(-Keep); 
                                TInvestorLikelihood(subj)= Keep;
                                TInvestorPlan(subj)=plan;
                                TInvestorToM(subj)=2*(level-1);                             
                                TInvestorAversion(subj) = risk;
                                TInvestorGuilt(subj) = guilt;
                                TInvestorTemp(subj) = tgrid(t);
                                TInvestorIrritability(subj) = irr;
                                TInvestorShift(subj) = shift;                                
                                end
                             end                                  
                         end
                     end 
                     for irr = 0:4
                   
                        if ne(level,3)
                         for guilt = 0:2
                             Keep =  fread(fid,1,'double');
                             su =  su+1; 
                             if exp(-Keep) > THold(subj)
                                THold(subj) = exp(-Keep);
                                MTrusteeLikelihood(subj)= Keep;
                                MTrusteePlan(subj)=plan;
                                MTrusteeToM(subj)=2*level-1;                             
                                MTrusteeAversion(subj) = risk;
                                MTrusteeGuilt(subj) = guilt;
                                MTrusteeTemp(subj) = tgrid(t);
                                MTrusteeIrritability(subj) = irr;
                                MTrusteeShift(subj) = shift;
                             end  
                             if ( shift == 0)
                                 ue = ue+1;
                                if exp(-Keep) > TirrHold(subj)                                 
                                TirrHold(subj) = exp(-Keep); 
                                ITrusteeLikelihood(subj)= Keep;
                                ITrusteePlan(subj)=plan;
                                ITrusteeToM(subj)=2*(level-1);                             
                                ITrusteeAversion(subj) = risk;
                                ITrusteeGuilt(subj) = guilt;
                                ITrusteeTemp(subj) = tgrid(t);
                                ITrusteeIrritability(subj) = irr;
                                ITrusteeShift(subj) = shift;                                
                                end
                             end                                
                             if (irr == 0 & shift == 0)
                                 u = u+1;
                                if exp(-Keep) > TavHold(subj)                                 
                                TavHold(subj) = exp(-Keep); 
                                AvTrusteeLikelihood(subj)= Keep;
                                AvTrusteePlan(subj)=plan;
                                AvTrusteeToM(subj)=2*(level-1);                             
                                AvTrusteeAversion(subj) = risk;
                                AvTrusteeGuilt(subj) = guilt;
                                AvTrusteeTemp(subj) = tgrid(t);
                                AvTrusteeIrritability(subj) = irr;
                                AvTrusteeShift(subj) = shift;                                
                                end
                             end      
                             if (irr == 0 & risk == 3 & t== 2 & shift == 0)
                                 o = o+1;
                                if exp(-Keep) > ToHold(subj)                                 
                                ToHold(subj) = exp(-Keep); 
                                OTrusteeLikelihood(subj)= Keep;
                                OTrusteePlan(subj)=plan;
                                OTrusteeToM(subj)=2*(level-1);                             
                                OTrusteeAversion(subj) = risk;
                                OTrusteeGuilt(subj) = guilt;
                                OTrusteeTemp(subj) = tgrid(t);
                                OTrusteeIrritability(subj) = irr;
                                OTrusteeShift(subj) = shift;                                
                                end
                             end        
                             if (irr == 0 & risk == 3 & shift == 0)
                                 he = he+1;
                                if exp(-Keep) > TtempHold(subj)                                 
                                TtempHold(subj) = exp(-Keep); 
                                TTrusteeLikelihood(subj)= Keep;
                                TTrusteePlan(subj)=plan;
                                TTrusteeToM(subj)=2*(level-1);                             
                                TTrusteeAversion(subj) = risk;
                                TTrusteeGuilt(subj) = guilt;
                                TTrusteeTemp(subj) = tgrid(t);
                                TTrusteeIrritability(subj) = irr;
                                TTrusteeShift(subj) = shift;                                
                                end
                             end                                
                         end      
                        end
                     end
                     
                 end


                     
                                
              end
          end
      end
   end
end
    fclose(fid);
end

X= (0);


n = 93;
subj = 0;
for j = 1:(length(X))
fid =fopen(['VarEstIrrZero.bin'],'r'); %SampleGenerative7
for i = 1:n
    %if i>1
    %    subj = subj+1;
    %end
    f = 0;
    s = 0;
    u =0;
    o=0;
    su = 0;
    w =0;
    g=0;    
    he = 0;
    ue =0;
    subj = subj+1;
   %Keep =  fread(fid,1,'double');
   for risk = 0:7
      % s = s+1;
      for t = 1:4
         % s = s+1;
          for plan = 1:4
              %s = s+1;
              for shift = 0:4
                 % s = s+1;
                 for level = 1:1
                     %s = s+1;
                     for irr = 0:4
                         %s = s+1;
                         for guilt = 0:2
                             s = s+1;
                             Keep =  fread(fid,1,'double');
                            % totalprob(subj) = totalprob(subj)+exp(-Keep);
                            if exp(-Keep) > IHold(subj)
                                IHold(subj) = exp(-Keep);
                                MInvestorLikelihood(subj)= Keep;
                                MInvestorPlan(subj)=plan;
                                MInvestorToM(subj)=2*(level-1);                             
                                MInvestorAversion(subj) = risk;
                                MInvestorGuilt(subj) = guilt;
                                MInvestorTemp(subj) = tgrid(t);
                                MInvestorIrritability(subj) = irr;
                                MInvestorShift(subj) = shift;
                            end
                             if ( shift == 0)
                                 ue = ue+1;
                                if exp(-Keep) > IirrHold(subj)                                 
                                IirrHold(subj) = exp(-Keep); 
                                IInvestorLikelihood(subj)= Keep;
                                IInvestorPlan(subj)=plan;
                                IInvestorToM(subj)=2*(level-1);                             
                                IInvestorAversion(subj) = risk;
                                IInvestorGuilt(subj) = guilt;
                                IInvestorTemp(subj) = tgrid(t);
                                IInvestorIrritability(subj) = irr;
                                IInvestorShift(subj) = shift;                                
                                end
                             end                              
                             if (irr == 0 & shift == 0)
                                 u = u+1;
                                if exp(-Keep) > IavHold(subj)                                 
                                IavHold(subj) = exp(-Keep); 
                                AvInvestorLikelihood(subj)= Keep;
                                AvInvestorPlan(subj)=plan;
                                AvInvestorToM(subj)=2*(level-1);                             
                                AvInvestorAversion(subj) = risk;
                                AvInvestorGuilt(subj) = guilt;
                                AvInvestorTemp(subj) = tgrid(t);
                                AvInvestorIrritability(subj) = irr;
                                AvInvestorShift(subj) = shift;                                
                                end
                             end      
                             if (irr == 0 & risk == 3 & t== 2 & shift == 0)
                                 o = o+1;
                                if exp(-Keep) > IoHold(subj)                                 
                                IoHold(subj) = exp(-Keep); 
                                OInvestorLikelihood(subj)= Keep;
                                OInvestorPlan(subj)=plan;
                                OInvestorToM(subj)=2*(level-1);                             
                                OInvestorAversion(subj) = risk;
                                OInvestorGuilt(subj) = guilt;
                                OInvestorTemp(subj) = tgrid(t);
                                OInvestorIrritability(subj) = irr;
                                OInvestorShift(subj) = shift;                                
                                end
                             end          
                             if (irr == 0 & risk == 3 & shift == 0)
                                 he = he+1;
                                if exp(-Keep) > ItempHold(subj)                                 
                                ItempHold(subj) = exp(-Keep); 
                                TInvestorLikelihood(subj)= Keep;
                                TInvestorPlan(subj)=plan;
                                TInvestorToM(subj)=2*(level-1);                             
                                TInvestorAversion(subj) = risk;
                                TInvestorGuilt(subj) = guilt;
                                TInvestorTemp(subj) = tgrid(t);
                                TInvestorIrritability(subj) = irr;
                                TInvestorShift(subj) = shift;                                
                                end
                             end                                
                         end
                        for guilt = 0:2
                             Keep =  fread(fid,1,'double');
                             su =  su+1; 
                             if exp(-Keep) > THold(subj)
                                THold(subj) = exp(-Keep);
                                MTrusteeLikelihood(subj)= Keep;
                                MTrusteePlan(subj)=plan;
                                MTrusteeToM(subj)=level-1;                             
                                MTrusteeAversion(subj) = risk;
                                MTrusteeGuilt(subj) = guilt;
                                MTrusteeTemp(subj) = tgrid(t);
                                MTrusteeIrritability(subj) = irr;
                                MTrusteeShift(subj) = shift;
                             end   
                             if ( shift == 0)
                                 ue = ue+1;
                                if exp(-Keep) > TirrHold(subj)                                 
                                TirrHold(subj) = exp(-Keep); 
                                ITrusteeLikelihood(subj)= Keep;
                                ITrusteePlan(subj)=plan;
                                ITrusteeToM(subj)=(level-1);                             
                                ITrusteeAversion(subj) = risk;
                                ITrusteeGuilt(subj) = guilt;
                                ITrusteeTemp(subj) = tgrid(t);
                                ITrusteeIrritability(subj) = irr;
                                ITrusteeShift(subj) = shift;                                
                                end
                             end                               
                             if (irr == 0 & shift == 0)
                                 u = u+1;
                                if exp(-Keep) > TavHold(subj)                                 
                                TavHold(subj) = exp(-Keep); 
                                AvTrusteeLikelihood(subj)= Keep;
                                AvTrusteePlan(subj)=plan;
                                AvTrusteeToM(subj)=(level-1);                             
                                AvTrusteeAversion(subj) = risk;
                                AvTrusteeGuilt(subj) = guilt;
                                AvTrusteeTemp(subj) = tgrid(t);
                                AvTrusteeIrritability(subj) = irr;
                                AvTrusteeShift(subj) = shift;                                
                                end
                             end      
                             if (irr == 0 & risk == 3 & t== 2 & shift == 0)
                                 o = o+1;
                                if exp(-Keep) > ToHold(subj)                                 
                                ToHold(subj) = exp(-Keep); 
                                OTrusteeLikelihood(subj)= Keep;
                                OTrusteePlan(subj)=plan;
                                OTrusteeToM(subj)=(level-1);                             
                                OTrusteeAversion(subj) = risk;
                                OTrusteeGuilt(subj) = guilt;
                                OTrusteeTemp(subj) = tgrid(t);
                                OTrusteeIrritability(subj) = irr;
                                OTrusteeShift(subj) = shift;                                
                                end
                             end   
                             if (irr == 0 & risk == 3 & shift == 0)
                                 he = he+1;
                                if exp(-Keep) > TtempHold(subj)                                 
                                TtempHold(subj) = exp(-Keep); 
                                TTrusteeLikelihood(subj)= Keep;
                                TTrusteePlan(subj)=plan;
                                TTrusteeToM(subj)=(level-1);                             
                                TTrusteeAversion(subj) = risk;
                                TTrusteeGuilt(subj) = guilt;
                                TTrusteeTemp(subj) = tgrid(t);
                                TTrusteeIrritability(subj) = irr;
                                TTrusteeShift(subj) = shift;                                
                                end
                             end                               
                         end      
                        
                     end
                     
                 end


                     
                                
              end
          end
      end
   end
end
    fclose(fid);
end




