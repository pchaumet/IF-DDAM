function plotimageincposreal(hlocal,event,ximage,imagem,imagexc,imageyc,imagezc)

val = get(hlocal,'Value');

switch val

case 1


figure(550)

set(550,'DefaultAxesFontName','Times')
set(550,'DefaultAxesFontSize',12)
set(550,'DefaultAxesFontWeight','Bold')
set(550,'DefaultTextfontName','Times')
set(550,'DefaultTextfontSize',12)
set(550,'DefaultTextfontWeight','Bold')
set(550,'Position',[0 0 1000 600])

subplot('position',[0.1 0.1 0.75 0.75])

imagesc(ximage,ximage,imagem.^2')
axis xy
caxis([min(min(imagem.^2)) max(max(imagem.^2))])
shading interp
axis equal
axis image
colorbar

xlabel('$x$','Interpreter','latex','Fontsize',18)
ylabel('$y$','Interpreter','latex','Fontsize',18)



case 2


figure(550)

set(550,'DefaultAxesFontName','Times')
set(550,'DefaultAxesFontSize',12)
set(550,'DefaultAxesFontWeight','Bold')
set(550,'DefaultTextfontName','Times')
set(550,'DefaultTextfontSize',12)
set(550,'DefaultTextfontWeight','Bold')
set(550,'Position',[0 0 1000 600])

subplot('position',[0.1 0.1 0.75 0.75])

imagesc(ximage,ximage,imagem')
axis xy
caxis([min(min(imagem)) max(max(imagem))])
shading interp
axis equal
axis image
colorbar

xlabel('$x$','Interpreter','latex','Fontsize',18)
ylabel('$y$','Interpreter','latex','Fontsize',18)


case 3

figure(550)

set(550,'DefaultAxesFontName','Times')
set(550,'DefaultAxesFontSize',12)
set(550,'DefaultAxesFontWeight','Bold')
set(550,'DefaultTextfontName','Times')
set(550,'DefaultTextfontSize',12)
set(550,'DefaultTextfontWeight','Bold')
set(550,'Position',[0 0 1000 600])

subplot('position',[0.1 0.1 0.8 0.8])
imagesc(ximage,ximage,abs(imagexc'))
axis xy

caxis([min(min((imagexc')))  max(max((imagexc')))])

shading interp
axis equal
axis image
colorbar
xlabel('$x$','Interpreter','latex','Fontsize',18)
ylabel('$y$','Interpreter','latex','Fontsize',18)

case 4

figure(550)

set(550,'DefaultAxesFontName','Times')
set(550,'DefaultAxesFontSize',12)
set(550,'DefaultAxesFontWeight','Bold')
set(550,'DefaultTextfontName','Times')
set(550,'DefaultTextfontSize',12)
set(550,'DefaultTextfontWeight','Bold')
set(550,'Position',[0 0 1000 600])


subplot('position',[0.1 0.1 0.8 0.8])
imagesc(ximage,ximage,abs(imageyc'))
axis xy
caxis([min(min(abs(imageyc'))) max(max(abs(imageyc')))])
  
shading interp
axis equal
axis image
colorbar
xlabel('$x$','Interpreter','latex','Fontsize',18)
ylabel('$y$','Interpreter','latex','Fontsize',18)

case 5

figure(550)

set(550,'DefaultAxesFontName','Times')
set(550,'DefaultAxesFontSize',12)
set(550,'DefaultAxesFontWeight','Bold')
set(550,'DefaultTextfontName','Times')
set(550,'DefaultTextfontSize',12)
set(550,'DefaultTextfontWeight','Bold')
set(550,'Position',[0 0 1000 600])

subplot('position',[0.1 0.1 0.8 0.8])

  imagesc(ximage,ximage,abs(imagezc'))

  axis xy
caxis([min(min(abs(imagezc'))) max(max(abs(imagezc')))])

shading interp
axis equal
axis image
colorbar
xlabel('$x$','Interpreter','latex','Fontsize',18)
ylabel('$y$','Interpreter','latex','Fontsize',18)

end;
