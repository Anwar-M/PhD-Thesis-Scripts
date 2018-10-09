function in = find_shadow(N_points,N_obs,tt,pt2,P_obs2)

%=====================================================================================================
%Find which points of a grid of observer are in the shadow of an object

%Input  - N_points , int	: number of points of the body
%		- N_obs , int		: number of receiver positions
%		- tt , int			: time step
%		- pt2, 2D array		: (N_points,3), coordinates of the points of the shielding body
%		- P_obs2 , 2D array	: (N_obs,3), coordinates of the receivers

%Output	- in, array			: 0 if in the shadow and 1 if in the light of the object for each receiver
 
%=====================================================================================================


	
for i=1:N_points 
	s = pt2(tt,i,1)^2+pt2(tt,i,2)^2;
	t = (1/(2*s))*( pt2(tt,i,3) + sqrt(pt2(tt,i,3)^2 + 4*s ) );
		
	pp(i,1) = t*pt2(tt,i,1);
	pp(i,2) = t*pt2(tt,i,2);
end

for i=1:N_obs
	s = P_obs2(tt,i,1)^2+P_obs2(tt,i,2)^2;
	t = (1/(2*s))*( P_obs2(tt,i,3) + sqrt(P_obs2(tt,i,3)^2 + 4*s ) );
		
	PP_obs(i,1) = t*P_obs2(tt,i,1);
	PP_obs(i,2) = t*P_obs2(tt,i,2);
end

xv = pp(:,1);
yv = pp(:,2);

xq =PP_obs(:,1);
yq =PP_obs(:,2);

in = inpolygon(xq,yq,xv,yv); %check which points are inside, value of csi

	
%Plot intersection
%figure(1)

%plot(xv,yv) % polygon
%axis equal

%hold on
%plot(xq(in),yq(in),'r+') % points inside
%plot(xq(~in),yq(~in),'bo') % points outside
%hold off