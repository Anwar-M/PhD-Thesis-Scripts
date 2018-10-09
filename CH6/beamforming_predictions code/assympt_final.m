% Diffraction integral method
% 09-06-2016
clear all;

tic 
z_pos = -1.5;
plane_observer = -7.5;
freq = 5000;%[50 63 80 100 125 160 200 250 315 400 500 630 800 1000 1250 1600 2000 2500 3150 4000 5000 6300 8000 10000];
n_freq = size(freq,2);
c = 344.0;
for ff=1:n_freq
	k_number(ff) = (2*pi*freq(ff))/c;
end
source_pos = [0 0 0];

%Read points from input file
in_file = fopen('coordinates_diff.txt','r');
if (in_file == -1) 
    error('File can''t be read');   
end
i=1;
while ~feof(in_file)
	
	line = fgets(in_file);
	[point_aux]=sscanf(line,'%f %f');
	pt(i,1) = point_aux(1)-source_pos(1); %Source position should be aligned with the axis origin
	pt(i,2) = point_aux(2)-source_pos(2);
	pt(i,3) = z_pos-source_pos(3);
	i=i+1; 
end
fclose(in_file);
N_points = size(pt,1);


%Generate grid of observers

%Number of observers in each direction
Nx =320;
Ny = 320;

%Limits of the grid
P1 = [-30 -30 0];
P2 = [-30  30 0];
P3 = [ 30  30 0];
P4 = [ 30 -30 0];

N_obs = Nx*Ny;

dx = (P3(1) - P1(1) )/(Nx-1);
dy = (P2(2) - P1(2) )/(Ny-1);

n=1;
for k=1:Ny
	for l=1:Nx
		P_obs(n,1:3) = [P1(1)+dx*(l-1) P1(2)+dy*(k-1) plane_observer];
		
		P_obs(n,1) = P_obs(n,1) - source_pos(1);
		P_obs(n,2) = P_obs(n,2) - source_pos(2);
		P_obs(n,3) = P_obs(n,3) - source_pos(3);
		
		n = n+1;
	end
end

clear P_obs;

P_obs(:,1)=0:0.01:6;
P_obs(:,2)=0;
P_obs(:,3)=-8;

N_obs=size(P_obs,1);

%Find shadow zone

for i=1:N_points
	s = pt(i,1)^2+pt(i,2)^2;
	t = (1/(2*s))*( pt(i,3) + sqrt(pt(i,3)^2 + 4*s ) );
	
	pp(i,1) = t*pt(i,1);
	pp(i,2) = t*pt(i,2);
end

for i=1:N_obs
	s = P_obs(i,1)^2+P_obs(i,2)^2;
	t = (1/(2*s))*( P_obs(i,3) + sqrt(P_obs(i,3)^2 + 4*s ) );
	
	PP_obs(i,1) = t*P_obs(i,1);
	PP_obs(i,2) = t*P_obs(i,2);
end

xv = pp(:,1);
yv = pp(:,2);

xq =PP_obs(:,1);
yq =PP_obs(:,2);



in = inpolygon(xq,yq,xv,yv); %check which points are inside, value of csi

%Exclude points on the edges
%for n=1:N_obs
%	for i=1:N_points
%		if(abs(pp(i,1)-PP_obs(n,1))<1e-4 && abs(pp(i,2)-PP_obs(n,2))<1e-4)
%			in(n) = 0;
%		end
%	end
%end
in(1)=1;


figure

plot(xv,yv) % polygon
axis equal

hold on
plot(xq(in),yq(in),'r+') % points inside
plot(xq(~in),yq(~in),'bo') % points outside
hold off
%-----------------
for ff=1:n_freq
	for n=1:N_obs

		for i=1:N_points-1
			
			L(i) = norm(pt(i+1,:)-pt(i,:));
		
			s(1) =0.0;
			s(i+1)=L(i)+s(i);
			
			e(i,:) = (pt(i+1,:)-pt(i,:))/L(i);
			y0(i,:) = pt(i,:) -s(i)*e(i,:);
			a(i,:) = y0(i,:);
			b(i,:) =y0(i,:) - P_obs(n,:);
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
				
				csi_s(i) = sqrt(k_number(ff)*abs( rho_s(i) + r_s(i) - norm(P_obs(n,:))));
				
				if(round(csi_s(i),3)==0)
					fresnel_s(i) = 1/2;
				else
					fresnel_s(i) = (1/(1 - exp(exp(-1i*pi/4)*2*sqrt(pi)*csi_s(i)))) + exp(1i*(csi_s(i)^2 + (pi/4)))/(2*sqrt(pi)*csi_s(i));
				end
			
				rho_r_s(i) = gamma(i) + (alpha(i)+beta(i))*s_s(i) + s_s(i)^2;
				

				f_s(i) =  A(i)/((r_s(i)*rho_s(i))*(rho_r_s(i) + rho_s(i)*r_s(i)));
				
				g(i) = (1/rho_s(i)) + (1/r_s(i)) - ((alpha(i)+s_s(i))^2)/(rho_s(i)^3) - ((beta(i) + s_s(i))^2)/(r_s(i))^3;
				
				h_s(i) = sqrt(k_number(ff)*(g(i)/2));
				
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
				
				csi_a(i) = sqrt(k_number(ff)*abs(rho_a(i) + ra(i) - r_s(i) - rho_s(i)));
				
				if(round(csi_a(i),4) == 0)
					fresnel_a(i) = 1/2;
				else
					fresnel_a(i) = (1/(1 - exp(exp(-1i*pi/4)*2*sqrt(pi)*csi_a(i)))) + exp(1i*(csi_a(i)^2 + (pi/4)))/(2*sqrt(pi)*csi_a(i));
				end
				
				rho_a_ra(i) = gamma(i) + (alpha(i) + beta(i))*sa(i) + (sa(i))^2;
			
				f_a(i) = A(i)/((ra(i)*rho_a(i))*(rho_a_ra(i)+rho_a(i)*ra(i)));
				
				ga(i) = (alpha(i) + sa(i) )/(rho_a(i)) + (beta(i)+sa(i))/(ra(i));
			
				ha(i) = k_number(ff)*ga(i)/(2*csi_a(i));
				
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
				
				csi_b(i) = sqrt(k_number(ff)*abs(rho_b(i) + rb(i) - r_s(i) - rho_s(i) ));
				
				if(round(csi_b(i),3) == 0)
					fresnel_b(i) = 1/2;
				else
					fresnel_b(i) = (1/(1 - exp(exp(-1i*pi/4)*2*sqrt(pi)*csi_b(i)))) + exp(1i*(csi_b(i)^2 + (pi/4)))/(2*sqrt(pi)*csi_b(i));
				end
				
				rho_b_rb(i) = gamma(i) + (alpha(i) + beta(i))*sb(i) + (sb(i))^2;
			
				f_b(i) = A(i)/((rb(i)*rho_b(i))*(rho_b_rb(i)+rho_b(i)*rb(i)));
				
				gb(i) = (alpha(i) + sb(i) )/(rho_b(i)) + (beta(i)+sb(i))/(rb(i));
			
				hb(i) = k_number(ff)*gb(i)/(2*csi_b(i));
				
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
			p_i = exp(1i*k_number(ff)*norm(P_obs(n,:)))/norm(P_obs(n,:));
		
		
			integral = sum(I_end + I_s);
		
		
			
			if( in(n)==0)
				integral2 = integral*exp(1i*k_number(ff)*norm(P_obs(n,:)));
			else
				integral2 = integral*exp(1i*k_number(ff)*norm(P_obs(n,:)))+p_i;
			end

			ps0=p_i-integral2;
			
			spl(n,ff) = 20*log10(abs((ps0)/p_i));
			fprintf('%d, spl %f\n',n,spl(n));
			%if(n==2866)
			%return;
			%end
			%if(isnan(spl(n,ff))==1)
			%	return;
			%end

	end

end
	toc
%Remove positive values outside the shadow region
%for ff=1:n_freq
%	for n=1:N_obs
%		if(spl(n,ff)>0)
%			spl(n,ff)=0.0;
%		end
%	end
%end

for n=1:N_obs
	OASPL(n) = -10*log10(10^(-spl(n,1)/10));
	for ff=2:n_freq
		
		if(spl(n,ff)<0)
			OASPL(n) = (10^(-spl(n,ff)/10) + 10^(-OASPL(n)/10));
		else
			OASPL(n) = -10*log10((1/n_freq)*OASPL(n));
		end
	end
end


toc
nn=1;
for(k=1:Ny)
	for l=1:Nx
		Attenuation(k,l) = OASPL(nn);
		nn = nn+1;
	end
end	
[X,Y] = meshgrid(P1(1):dx:P4(1),P1(2):dy:P2(2));

colormap(jet(20));
pcolor(X,Y,Attenuation);
xlabel('x(m)');   
ylabel('y(m)');
set(gca,'color','none')
shading interp;

bar=colorbar;
%V=[-15 0];      %Define contour values
%caxis(V);   
 ptt=[pt]; 
hold on
plot(ptt(:,1)+source_pos(1),ptt(:,2)+source_pos(2),'w-','LineWidth',2); %Define the shielding object boundary
plot(source_pos(1),source_pos(2),'m.','markersize',15);
axis equal; 
hold off

