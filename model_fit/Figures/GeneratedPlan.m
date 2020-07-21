X = [0];% 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23];

IDisapp = zeros(1,4);
n = 93;
for j = 1:(length(X))
fid =fopen(['ExploitGen' '.bin'],'r'); %SampleGenerative7 21Output Averse7
for i = 1:n
 for t = 1:10    
    Investor(i+(j-1)*n,t)=fread(fid,1,'int32');
    Trustee(i+(j-1)*n,t)=fread(fid,1,'int32');
 end
end
    fclose(fid);
end

for j = 1:(length(X))
fid =fopen(['ImpulsiveGen' '.bin'],'r'); %SampleGenerative7 21Output Averse7
for i = 1:n
 for t = 1:10    
    XInvestor(i+(j-1)*n,t)=fread(fid,1,'int32');
    XTrustee(i+(j-1)*n,t)=fread(fid,1,'int32');
 end
end
    fclose(fid);
end

for j = 1:(length(X))
fid =fopen(['PlanMismatchGen' '.bin'],'r'); %SampleGenerative7 21Output Averse7
for i = 1:n
 for t = 1:10    
    YInvestor(i+(j-1)*n,t)=fread(fid,1,'int32');
    YTrustee(i+(j-1)*n,t)=fread(fid,1,'int32');
 end
end
    fclose(fid);
end