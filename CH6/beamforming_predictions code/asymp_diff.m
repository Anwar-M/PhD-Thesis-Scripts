function [I_s,I_end] = asymp_diff(tt,pt2,N_points,P_obs2,n,k_number)

%==================================================================================
%Applies the assymptotic expansion of the diffraction integral method

%Input  - tt , int			: time step
%		- pt2, 2D array		: (N_points,3) coordinates of the points of the body 
%		- N_poins , int		: number of points in the body
%		- P_obs2 , 2D array	: (N observers,3) coordinates of the observers
%		- n , int			: observer considered
%		- k_number, array	: wave number for each frequency

%Output	- I_s, img			: contribution of the stationary phase
%		- I_end, img		: contribution of the end points
%==================================================================================

for i=1:N_points-1
	
	L(i) = sqrt( (pt2(tt,i+1,1)-pt2(tt,i,1))^2 + (pt2(tt,i+1,2)-pt2(tt,i,2))^2 +(pt2(tt,i+1,3)-pt2(tt,i,3))^2);
	
	s(1) =0.0;
	s(i+1)=L(i)+s(i);
	
	e(i,:) = (pt2(tt,i+1,:)-pt2(tt,i,:))/L(i);
				
	y0(i,1) = pt2(tt,i,1) - s(i)*e(i,1);
	y0(i,2) = pt2(tt,i,2) - s(i)*e(i,2);
	y0(i,3) = pt2(tt,i,3) - s(i)*e(i,3);
	
	a(i,:) = y0(i,:);
				
	b(i,1) =y0(i,1) - P_obs2(tt,n,1);
	b(i,2) =y0(i,2) - P_obs2(tt,n,2);
	b(i,3) =y0(i,3) - P_obs2(tt,n,3);
	
	a2(i) = dot(a(i,:),a(i,:));
	b2(i) = dot(b(i,:),b(i,:));
	alpha(i) = dot(a(i,:),e(i,:));
	beta(i) = dot(b(i,:),e(i,:));
	gamma(i) = dot(a(i,:),b(i,:));
			
	u(i,:) = cross(a(i,:),b(i,:));
	aux_v = a(i,:)-b(i,:);
	v(i,:) = cross(e(i,:),aux_v);
	
	u2(i) = dot(u(i,:),u(i,:));
	v2(i) = dot(v(i,:),v(i,:));
			
	uv(i) = dot(u(i,:),v(i,:));
	ue(i) = dot(u(i,:),e(i,:));
				
	w(i,:) = cross(a(i,:),e(i,:));
	z(i,:) = v(i,:) + w(i,:);
				
	delta_s(i) = dot(w(i,:),w(i,:));
	delta_p(i) = dot(z(i,:),z(i,:));
	
end


for i=1:N_points-1

	% Calculate stationary phase contribution
					
	A(i) = ue(i)/(4*pi);
	sa(i) = s(i);
	sb(i) = s(i+1);
			
	delta_ss(i) = sqrt(delta_s(i));
	delta_pp(i) = sqrt(delta_p(i));
	
	s_s(i) = -(delta_ss(i)*beta(i) + delta_pp(i)*alpha(i))/( delta_ss(i) + delta_pp(i) );
			
	rho_s(i) = sqrt(a2(i) + 2*alpha(i)*s_s(i) + s_s(i)^2);
					
	r_s(i) = sqrt(b2(i) + 2*beta(i)*s_s(i) + s_s(i)^2);
					
	csi_s(i) = sqrt(k_number*abs( rho_s(i) + r_s(i) - sqrt(P_obs2(tt,n,1)^2 + P_obs2(tt,n,2)^2 + P_obs2(tt,n,3)^2)));
					
	if(round(csi_s(i),3)==0)
		fresnel_s(i) = 1/2;
	else
		fresnel_s(i) = (1/(1 - exp(exp(-1i*pi/4)*2*sqrt(pi)*csi_s(i)))) + exp(1i*(csi_s(i)^2 + (pi/4)))/(2*sqrt(pi)*csi_s(i));
	end
			
	rho_r_s(i) = gamma(i) + (alpha(i)+beta(i))*s_s(i) + s_s(i)^2;
	
	f_s(i) =  A(i)/((r_s(i)*rho_s(i))*(rho_r_s(i) + rho_s(i)*r_s(i)));
					
	g(i) = (1/rho_s(i)) + (1/r_s(i)) - ((alpha(i)+s_s(i))^2)/(rho_s(i)^3) - ((beta(i) + s_s(i))^2)/(r_s(i))^3;
					
	h_s(i) = sqrt(k_number*(g(i)/2));
					
	G_s(i)=f_s(i)/h_s(i);
				
	if(isnan(G_s(i)))
		G_s(i)=0;
	end
				
	I_s(i) = 2*pi*csi_s(i)*fresnel_s(i)*G_s(i);
					
	indicator = (s_s(i)>sa(i))*(s_s(i)<sb(i));
					
	I_s(i)=indicator*I_s(i);
	
	%Calculate starting point sa, contribution
					
	rho_a(i) = sqrt(a2(i) + 2*alpha(i)*sa(i) + sa(i)*sa(i) );
					
	ra(i) =  sqrt(b2(i) + 2*beta(i)*sa(i) + sa(i)*sa(i) );
					
	csi_a(i) = sqrt(k_number*abs(rho_a(i) + ra(i) - r_s(i) - rho_s(i)));
					
	if(round(csi_a(i),4) == 0)
		fresnel_a(i) = 1/2;
	else
		fresnel_a(i) = (1/(1 - exp(exp(-1i*pi/4)*2*sqrt(pi)*csi_a(i)))) + exp(1i*(csi_a(i)^2 + (pi/4)))/(2*sqrt(pi)*csi_a(i));
	end

	rho_a_ra(i) = gamma(i) + (alpha(i) + beta(i))*sa(i) + (sa(i))^2;
			
	f_a(i) = A(i)/((ra(i)*rho_a(i))*(rho_a_ra(i)+rho_a(i)*ra(i)));
					
	ga(i) = (alpha(i) + sa(i) )/(rho_a(i)) + (beta(i)+sa(i))/(ra(i));
				
	ha(i) = k_number*ga(i)/(2*csi_a(i));
					
	if(round(s_s(i)/sa(i),3)==1)
		ind=1;
	else
		ind=0;
	end
					
	G_a1(i) = (f_s(i)/h_s(i))*(sa(i)==s_s(i));

	G_a2(i) = (f_a(i)/ha(i))*(1-(sa(i)==s_s(i)));
					
	if(isnan(G_a1(i)))
		G_a1(i)=0.0;
	end
	if(isnan(G_a2(i)))
		G_a2(i)=0.0;
	end
	G_a(i)=G_a1(i)+G_a2(i);
				
			
	Ia(i) = fresnel_a(i)*G_a(i);


	%Calculate starting point sb, contribution
					
	rho_b(i) = sqrt(a2(i) + 2*alpha(i)*sb(i) + sb(i)*sb(i) );
			
	rb(i) =  sqrt(b2(i) + 2*beta(i)*sb(i) + sb(i)*sb(i) );
					
	csi_b(i) = sqrt(k_number*abs(rho_b(i) + rb(i) - r_s(i) - rho_s(i) ));
					
	if(round(csi_b(i),3) == 0)
		fresnel_b(i) = 1/2;
	else
		fresnel_b(i) = (1/(1 - exp(exp(-1i*pi/4)*2*sqrt(pi)*csi_b(i)))) + exp(1i*(csi_b(i)^2 + (pi/4)))/(2*sqrt(pi)*csi_b(i));
	end
				
	rho_b_rb(i) = gamma(i) + (alpha(i) + beta(i))*sb(i) + (sb(i))^2;
			
	f_b(i) = A(i)/((rb(i)*rho_b(i))*(rho_b_rb(i)+rho_b(i)*rb(i)));
					
	gb(i) = (alpha(i) + sb(i) )/(rho_b(i)) + (beta(i)+sb(i))/(rb(i));
				
	hb(i) = k_number*gb(i)/(2*csi_b(i));
					
	if(round(s_s(i)/sb(i),4)==1)
		ind=1;
	else
		ind=0;
	end
				
	G_b1(i) = (f_s(i)/h_s(i))*(sb(i)==s_s(i));
				
	G_b2(i) = (f_b(i)/hb(i))*(1-(sb(i)==s_s(i)));
				
	if(isnan(G_b1(i)))
		G_b1(i)=0.0;
	end
	if(isnan(G_b2(i)))
		G_b2(i)=0.0;
	end
				
	G_b(i)=G_b1(i)+G_b2(i);
	Ib(i) = fresnel_b(i)*G_b(i);
					
	I_end(i) = 2*pi*abs(csi_s(i))*fresnel_s(i)*(Ia(i)-Ib(i));

	
end