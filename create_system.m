function create_system(phi)
    addpath '../distmesh'
    addpath '../'

    height_shift = 0;
    play = 0;
    
    %% Cell Body
    cell_radius = 2.5;
    eps = 0.8;
    if 1
        fd=@(p) dsphere(p,0,0,0,cell_radius);
        [mp,mt]=distmeshsurface(fd,@huniform,eps,1.1*[-3,-3,-3;3,3,3],10000);
        save sphere mp mt
    else
        [X,Y,Z] = sphere(2*round(2*cell_radius/eps));
        X = cell_radius*X;
        Y = cell_radius*Y;
        Z = cell_radius*Z;
        a = surf2patch(X,Y,Z,'triangles');
        mt = a.faces;
        mp = a.vertices;
    end
    mv = 0*mt;

    %% Flagellum
    flagellum_radius = 0.25;
    length = 6;
    width = 2.5;
    k = 2*pi;
    
    s = linspace(0,1,100)';
    sx = [length*s,width*(1-exp(-k*s).^2).*cos(k*s-phi),...
                width*(1-exp(-k*s).^2).*sin(k*s-phi)];
    sv = [0*s,width*(1-exp(-k*s).^2).*sin(k*s-phi),...
                -width*(1-exp(-k*s).^2).*cos(k*s-phi)];
    sx(:,1) = sx(:,1)+cell_radius+0.1;
    sw = flagellum_radius*ones(size(sx,1),1);
    
    if ~play
        %% Remove duplicates
        for i=size(sx,1):-1:1
            for j=1:i-1
                if sum(abs(sx(i,:)-sx(j,:)))<0.001
                    sx(i,:) = [];
                    sv(i,:) = [];
                    sw(i,:) = [];
                end
            end
        end
    end

    %% Shift
    mp(:,3) = mp(:,3)+height_shift;
    sx(:,3) = sx(:,3)+height_shift;
    
    %% Plot
    clf()
    patch('faces',mt,'vertices',mp,'facecol',cont_colour(0.5),'edgecol','k','facealpha',1)
    hold on
    plot3(sx(:,1),sx(:,2),sx(:,3),'x')
    hold on
    quiver3(sx(:,1),sx(:,2),sx(:,3),sv(:,1),sv(:,2),sv(:,3),0)
    view(45,20)
    axis equal
    % plot3(evu_x,evu_y,evu_z,'.')
    view(0,0)
    drawnow()
    if play
        disp('play=1 : not creating files!')
        return
    end
 
    
    %% Evaulate velocity at
    coarse_dx = 1.5;
    fine_dx = 0.35;
    
    xx = -10:coarse_dx:25;
    yy = -15:coarse_dx:15;
    zz = horzcat(-wrev(coarse_dx:coarse_dx:10),0:coarse_dx:10)+height_shift;
    
    
    [xx,yy,zz] = meshgrid(xx,yy,zz);
    evu_x = reshape(xx,1,numel(xx));
    evu_y = reshape(yy,1,numel(yy));
    evu_z = reshape(zz,1,numel(zz));
    
    xx = 0:fine_dx:12.5;
    yy = -7.5:fine_dx:7.5;
    zz = horzcat(-wrev(fine_dx:fine_dx:1.5),0:fine_dx:1.5)+height_shift;
    [xx,yy,zz] = meshgrid(xx,yy,zz);
   
    evu_x(end+1:end+numel(xx)) = reshape(xx,1,numel(xx));
    evu_y(end+1:end+numel(xx)) = reshape(yy,1,numel(yy));
    evu_z(end+1:end+numel(xx)) = reshape(zz,1,numel(zz));
    for i=numel(evu_x):-1:1 % Remove duplicates
        for j=1:(i-1)
            if (evu_x(i)==evu_x(j)) && ...
                    (evu_y(i)==evu_y(j)) && ...
                    (evu_z(i)==evu_z(j))
                evu_x(i) = [];
                evu_y(i) = [];
                evu_z(i) = [];
            end
        end
    end
    
    %% Output to txt
    % Output Mesh
    disp('Outputting mesh.')
    fid = fopen( 'system.txt', 'wt' );
    fprintf(fid, '%i\n',size(mt,1));
    for i=1:size(mt,1)
        fprintf(fid, '%i %i %i\n',mt(i,1),mt(i,2),mt(i,3));
    end
    fprintf(fid, '%i\n',size(mp,1));
    for i=1:size(mp,1)
        fprintf(fid, '%f %f %f\n', mp(i,1),mp(i,2),mp(i,3));
    end
    fprintf(fid, '%i\n',size(mv,1));
    for i=1:size(mv,1)
        fprintf(fid, '%f %f %f\n', mv(i,1),mv(i,2),mv(i,3));
    end
    % Output Splines
    fprintf(fid, '%i\n',size(sx,1));
    for i=1:size(sx,1)
        fprintf(fid, '%i %i %i\n',sx(i,1),sx(i,2),sx(i,3));
    end
    fprintf(fid, '%i\n',size(sv,1));
    for i=1:size(sv,1)
        fprintf(fid, '%f %f %f\n', sv(i,1),sv(i,2),sv(i,3));
    end
    fprintf(fid, '%i\n',size(sw,1));
    for i=1:size(sw,1)
        fprintf(fid, '%f\n', sw(i));
    end
    fprintf(fid, '%i\n',size(evu_x,2));
    for i=1:size(evu_x,2)
        fprintf(fid, '%f %f %f\n', evu_x(i),evu_y(i),evu_z(i));
    end
    fclose(fid);
    save system
end