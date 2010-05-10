function outstring=protect(instring)
% protects eg. greek characters in a string
% whereever there is a backslash, it adds another

j=0;
for i=1:length(instring)
  j=j+1;
  outstring(j)=instring(i);
  if (instring(i)=='\') 
    j=j+1; 
    outstring(j)='\'; 
  end
end