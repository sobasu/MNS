 %% Output is "connectivity"
 %% "connectivity" lists [element number on Gamma1    nonzero dof on Omg2    value at first quadrature point    value at second quadrature point]
 
 load('mortar_deep_4_patches_curved_data')
  
  dyadic_power=3;
  k2=-1; 
  k1=0;
  n_ref_srf = 2^dyadic_power + k2;
  degree_of_spline_basis_srf=3;
  regularity_of_spline_basis_srf=2;
  
%% Creating curvilinear mesh for srf
 
    knts=[zeros(1, degree_of_spline_basis_srf) linspace(0,1, n_ref_srf) ones(1, degree_of_spline_basis_srf)];
    
    p=degree_of_spline_basis_srf;
    
    ncp = numel(knts) - (degree_of_spline_basis_srf+1);
  
  Amplitude=2^(-dyadic_power-2);
  
  t=linspace(0,1, ncp); %n_ref_srf + 1);
  yt1=sin(pi*t)*Amplitude;
   yt2=sin(2*pi*t)*Amplitude;
   
  size_of_inner_square=1/2;
  scaling_factor=size_of_inner_square + 2 *epsilon_by_h/ (2^dyadic_power + k1 ); %2/sqrt(40); 
  
  ybnd1 = linspace(1/2*(1-scaling_factor), 1/2*(1 + scaling_factor), ncp); %n_ref_srf + 1);
  xbnd1 = 1/2*(1-scaling_factor) + yt1;
  ybnd2 = ybnd1;
  xbnd2 = 1/2*(1 + scaling_factor) + yt2;
  
  xbnd3 = ybnd1;
  ybnd3 = 1/2*(1-scaling_factor) + yt2;
  xbnd4 = ybnd1;
  ybnd4 = 1/2*(1+scaling_factor) + yt1;
  
  zbnd=zeros(size(t));
  
  crv1 = nrbmak([xbnd1; ybnd1], knts);
  crv2 = nrbmak([xbnd2; ybnd2], knts);
  crv3 = nrbmak([xbnd3; ybnd3], knts);
  crv4 = nrbmak([xbnd4; ybnd4], knts);

  srf = nrbcoons(crv1, crv2, crv3, crv4);
  
  problem_data.geo_name=srf; 
    problem_data.c_diff  = @(x, y) ones(size(x));

%     dyadic_srf=dyadic_degree;
    method_data_srf.degree=repmat(degree_of_spline_basis_srf,1,2);
    method_data_srf.regularity = repmat(regularity_of_spline_basis_srf,1,2);
    method_data_srf.nquad      = [method_data_srf.degree+1];       % Points for the Gaussian quadrature rule
    method_data_srf.nsub       = [2 2].^0;       % Number of subdivisions
    

    [geometry_srf, msh_srf, space_srf] = ...
              get_mesh_space (problem_data, method_data_srf) ;
          
    subplot(1,2,1)
    nrbplot(srf,msh_srf.nel_dir);
    view(0,90); shg
    
    subplot(1,2,2)
    [X,Y]=meshgrid(unique(knts), unique(knts))
    mesh(X,Y, zeros(size(X)) )
    view(0,90); axis square; shg
%%

colorside={'b','r','k','m'};


plotsmbol='o .';
cnn=[];
connectivity=[];

for iside=[1 2 3 4]
    iparit=1;
    
qnn=cell2mat(...
    msh_Omg1.msh_patch{ edge_map.patch( find(edge_map.dOmg2==iside) )}. ...
    boundary(edge_map.dOmg1_inner_local(find(edge_map.dOmg2==iside))).qn);

    for inel=1 : size(qnn, 2)
       for iqn=1 : size(qnn, 1)
           if           edge_map.dOmg1_inner_local(find(edge_map.dOmg2==iside))==1 % Find the local numering of a patch side corresponding to side of srf
                pl = nrbeval(nrb(edge_map.patch( find(edge_map.dOmg2==iside) )), [0 qnn(iqn, inel)]); %Evaluate physical cordinates of each Omg1 boundary quadrature point
                loc{iside,inel,iqn} = flipud(nrbinverse(srf, pl)); % Evluate in srf parametric coordinates each Omg1 boundary quadrature point 
                fndspn{iside,inel,iqn, 1}= findspan(numel(knts)-degree_of_spline_basis_srf-1,...
                    degree_of_spline_basis_srf, loc{iside,inel,iqn}(1) ,knts);
                tval{iside,inel,iqn, 1} = basisfunder (fndspn{iside,inel,iqn, 1}, degree_of_spline_basis_srf, loc{iside,inel,iqn}(1) ,knts, 0);
                fndspn{iside,inel,iqn, 2}= findspan(numel(knts)-degree_of_spline_basis_srf-1,...
                    degree_of_spline_basis_srf, loc{iside,inel,iqn}(2) ,knts);
                tval{iside,inel,iqn, 2} = basisfunder (fndspn{iside,inel,iqn, 2}, degree_of_spline_basis_srf, loc{iside,inel,iqn}(2) ,knts, 0);
                    cbf{iside,inel,iqn, 1} = numbasisfun (fndspn{iside,inel,iqn, 1}, loc{iside,inel,iqn}(1), p, knts) + 1;
                        cbf{iside,inel,iqn, 2} = numbasisfun (fndspn{iside,inel,iqn, 2}, loc{iside,inel,iqn}(2), p, knts) + 1;
                 %two_dim_span=       
            elseif            edge_map.dOmg1_inner_local(find(edge_map.dOmg2==iside))==2
                pl = nrbeval(nrb(edge_map.patch( find(edge_map.dOmg2==iside) )), [1 qnn(iqn, inel)]);
                loc{iside,inel,iqn} = flipud(nrbinverse(srf, pl));
                fndspn{iside,inel,iqn, 1}= findspan(numel(knts)-degree_of_spline_basis_srf-1,...
                    degree_of_spline_basis_srf, loc{iside,inel,iqn}(1) ,knts);
                tval{iside,inel,iqn, 1} = basisfunder (fndspn{iside,inel,iqn, 1}, degree_of_spline_basis_srf, loc{iside,inel,iqn}(1) ,knts, 0);
                fndspn{iside,inel,iqn, 2}= findspan(numel(knts)-degree_of_spline_basis_srf-1,...
                    degree_of_spline_basis_srf, loc{iside,inel,iqn}(2) ,knts);
                tval{iside,inel,iqn, 2} = basisfunder (fndspn{iside,inel,iqn, 2}, degree_of_spline_basis_srf, loc{iside,inel,iqn}(2) ,knts, 0);
                    cbf{iside,inel,iqn, 1} = numbasisfun (fndspn{iside,inel,iqn, 1}, loc{iside,inel,iqn}(1), p, knts) + 1;
                        cbf{iside,inel,iqn, 2} = numbasisfun (fndspn{iside,inel,iqn, 2}, loc{iside,inel,iqn}(2), p, knts) + 1;
            elseif            edge_map.dOmg1_inner_local(find(edge_map.dOmg2==iside))==3
                pl = nrbeval(nrb(edge_map.patch( find(edge_map.dOmg2==iside) )), [qnn(iqn, inel) 0]);
                loc{iside,inel,iqn} = flipud(nrbinverse(srf, pl));
                fndspn{iside,inel,iqn, 1}= findspan(numel(knts)-degree_of_spline_basis_srf-1,...
                    degree_of_spline_basis_srf, loc{iside,inel,iqn}(1) ,knts);
                tval{iside,inel,iqn, 1} = basisfunder (fndspn{iside,inel,iqn, 1}, degree_of_spline_basis_srf, loc{iside,inel,iqn}(1) ,knts, 0);
                fndspn{iside,inel,iqn, 2}= findspan(numel(knts)-degree_of_spline_basis_srf-1,...
                    degree_of_spline_basis_srf, loc{iside,inel,iqn}(2) ,knts);
                tval{iside,inel,iqn, 2} = basisfunder (fndspn{iside,inel,iqn, 2}, degree_of_spline_basis_srf, loc{iside,inel,iqn}(2) ,knts, 0);
                    cbf{iside,inel,iqn, 1} = numbasisfun (fndspn{iside,inel,iqn, 1}, loc{iside,inel,iqn}(1), p, knts) + 1;
                        cbf{iside,inel,iqn, 2} = numbasisfun (fndspn{iside,inel,iqn, 2}, loc{iside,inel,iqn}(2), p, knts) + 1;
            else
                pl = nrbeval(nrb(edge_map.patch( find(edge_map.dOmg2==iside) )), [qnn(iqn, inel) 1]);
                loc{iside,inel,iqn}=flipud(nrbinverse(srf, pl));
                fndspn{iside,inel,iqn, 1}= findspan(numel(knts)-degree_of_spline_basis_srf-1,...
                    degree_of_spline_basis_srf, loc{iside,inel,iqn}(1) ,knts);
                tval{iside,inel,iqn, 1} = basisfunder (fndspn{iside,inel,iqn, 1}, degree_of_spline_basis_srf, loc{iside,inel,iqn}(1) ,knts, 0);
                fndspn{iside,inel,iqn, 2}= findspan(numel(knts)-degree_of_spline_basis_srf-1,...
                    degree_of_spline_basis_srf, loc{iside,inel,iqn}(2) ,knts);
                tval{iside,inel,iqn, 2} = basisfunder (fndspn{iside,inel,iqn, 2}, degree_of_spline_basis_srf, loc{iside,inel,iqn}(2) ,knts, 0);
                    cbf{iside,inel,iqn, 1} = numbasisfun (fndspn{iside,inel,iqn, 1}, loc{iside,inel,iqn}(1), p, knts) + 1;
                        cbf{iside,inel,iqn, 2} = numbasisfun (fndspn{iside,inel,iqn, 2}, loc{iside,inel,iqn}(2), p, knts) + 1;
           iparit=-iparit;
           end
           
       clear t1 t2 tc1 tc2

       for it=1:size(tval{iside,inel,iqn, 1}, 3)
           t1(it)=tval{iside,inel,iqn, 1}(it);
           t2(it)=tval{iside,inel,iqn, 2}(it);
       end
       tval_2d{iside,inel,iqn}=t1(:)*t2(:)';
       
       for it1=1:size(tval{iside,inel,iqn, 1}, 3);
           for it2=1:size(tval{iside,inel,iqn, 1}, 3);
               tcbf(it1, it2) = sub2ind(space_srf.ndof_dir, cbf{iside,inel,iqn, 1}(it1), cbf{iside,inel,iqn, 2}(it2));
           end
       end
       tcbf_2d{iside,inel,iqn}=tcbf;
        
       subplot(1,2,1)
       hold all
            plot3(pl(1), pl(2), pl(3), strcat(colorside{iside}, plotsmbol(iparit+2)))
       hold off
       
       subplot(1,2,2)
       hold all
            plot3(loc{iside,inel,iqn}(1), loc{iside,inel,iqn}(2), 0, strcat(colorside{iside}, plotsmbol(iparit+2)))
       hold off
       
       cnn=[cnn; iside inel iqn NaN tcbf(:)'];
       
     end % enf for loop over iqn
       
       U=[];
       for iqn=1:size(qnn, 1)
           U=union(U, tcbf_2d{iside, inel,iqn}(:));
       end
       
       conn=zeros(numel(U), size(qnn, 1));
       for iu=1:numel(U)
           for iqn=1:size(qnn, 1)
                try
                    conn(iu, iqn) =  tval_2d{iside,inel,iqn}(find(tcbf_2d{iside,inel,iqn}==U(iu)));
                end
           end
       end
       
        conn=[repmat(inel,[numel(U(:)),1]) U(:) conn];
        connectivity=[connectivity; conn];
    end
     
end

