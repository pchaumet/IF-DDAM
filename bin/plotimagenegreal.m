function plotimagenegreal(hlocal,event,ximage,imagem,imagexc,imageyc,imagezc,nprint)

val = get(hlocal,'Value');

switch val

case 1


figure(700)

set(700,'DefaultAxesFontName','Times')
set(700,'DefaultAxesFontSize',12)
set(700,'DefaultAxesFontWeight','Bold')
set(700,'DefaultTextfontName','Times')
set(700,'DefaultTextfontSize',12)
set(700,'DefaultTextfontWeight','Bold')
set(700,'Position',[0 0 1000 600])

subplot('position',[0.1 0.1 0.8 0.8])

imagesc(ximage,ximage,imagem.^2')
axis xy
caxis([ 0 max(max(imagem.^2))])
shading interp
axis equal
axis image
colorbar

xlabel('$x$','Interpreter','latex','Fontsize',18)
ylabel('$y$','Interpreter','latex','Fontsize',18)


case 2


figure(700)

set(700,'DefaultAxesFontName','Times')
set(700,'DefaultAxesFontSize',12)
set(700,'DefaultAxesFontWeight','Bold')
set(700,'DefaultTextfontName','Times')
set(700,'DefaultTextfontSize',12)
set(700,'DefaultTextfontWeight','Bold')
set(700,'Position',[0 0 1000 600])

subplot('position',[0.1 0.1 0.8 0.8])

imagesc(ximage,ximage,imagem')
axis xy
caxis([ 0 max(max(imagem))])
shading interp
axis equal
axis image
colorbar

xlabel('$x$','Interpreter','latex','Fontsize',18)
ylabel('$y$','Interpreter','latex','Fontsize',18)


case 3

figure(700)

set(700,'DefaultAxesFontName','Times')
set(700,'DefaultAxesFontSize',12)
set(700,'DefaultAxesFontWeight','Bold')
set(700,'DefaultTextfontName','Times')
set(700,'DefaultTextfontSize',12)
set(700,'DefaultTextfontWeight','Bold')
set(700,'Position',[0 0 1000 600])

subplot('position',[0.1 0.1 0.8 0.8])

imagesc(ximage,ximage,abs(imagexc'))
axis xy

caxis([ 0 max(max(abs(imagexc')))])

shading interp
axis equal
axis image
colorbar
xlabel('$x$','Interpreter','latex','Fontsize',18)
ylabel('$y$','Interpreter','latex','Fontsize',18)

case 4

figure(700)

set(700,'DefaultAxesFontName','Times')
set(700,'DefaultAxesFontSize',12)
set(700,'DefaultAxesFontWeight','Bold')
set(700,'DefaultTextfontName','Times')
set(700,'DefaultTextfontSize',12)
set(700,'DefaultTextfontWeight','Bold')
set(700,'Position',[0 0 1000 600])

  
subplot('position',[0.1 0.1 0.8 0.8])

imagesc(ximage,ximage,abs(imageyc'))
axis xy
caxis([ 0 max(max(abs(imageyc')))])
  
shading interp
axis equal
axis image
colorbar
xlabel('$x$','Interpreter','latex','Fontsize',18)
ylabel('$y$','Interpreter','latex','Fontsize',18)

case 5

figure(700)

set(700,'DefaultAxesFontName','Times')
set(700,'DefaultAxesFontSize',12)
set(700,'DefaultAxesFontWeight','Bold')
set(700,'DefaultTextfontName','Times')
set(700,'DefaultTextfontSize',12)
set(700,'DefaultTextfontWeight','Bold')
set(700,'Position',[0 0 1000 600])


subplot('position',[0.1 0.1 0.8 0.8])
imagesc(ximage,ximage,abs(imagezc'))

axis xy
caxis([ 0 max(max(abs(imagezc')))])

shading interp
axis equal
axis image
colorbar
xlabel('$x$','Interpreter','latex','Fontsize',18)
ylabel('$y$','Interpreter','latex','Fontsize',18)

end;

if (nprint == 1)
print('-f700','imagenegwf','-depsc')
end
