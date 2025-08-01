
%% THIS FUNCTION RETURNS U,V VALUES AT T-CARRIER PTS

function [Uint,Vint] = fun_int(U,V)


% Get x,y dimensions of U,V

nx = size(U,1);
ny = size(U,2);


% Interpolate/extrapolate U values at T-carrier pts 

for i=1:nx-1
for j=1:ny
    Uint(i,j,:) = (U(i,j,:)+U(i+1,j,:))/2;
end
end

for j=1:ny
    Uint(nx,j,:) = U(nx,j,:);
end


% Interpolate/extrapolate V values at T-carrier pts

for i=1:nx
for j=1:ny-1
    Vint(i,j,:) = (V(i,j,:)+V(i,j+1,:))/2;
end
end

for i=1:nx
    Vint(i,ny,:) = V(i,ny,:);
end