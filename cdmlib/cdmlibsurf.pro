######################################################################
# Automatically generated by qmake (2.01a) Mon Dec 27 12:07:07 2010
######################################################################

TEMPLATE 	= 	lib

VERSION         =       0.1.7

TARGET 		=       cdmlibsurf

CONFIG          +=      staticlib warn_on

DEPENDPATH 	+= .

DESTDIR      	= lib

QT      	+=

DEFINES 	+=      CDMVERSION=\\\"$$VERSION\\\"

DEFINES 	+= 	QT_NO_DEBUG_OUTPUT

QMAKE_CC        =       gfortran -Warray-bounds -fcray-pointer -w

QMAKE_CFLAGS    += -fopenmp -lfftw3_omp -lfftw3 -lm 

QMAKE_CFLAGS_RELEASE    = -O3 -w

QMAKE_CFLAGS_THREAD =

HEADERS 	+=      cdmlibsurf.h

SOURCES		+= aleatoire.f \
                bessel.f \
                bicgstarplus.f \
		cdmlibsurf.f \
		champcircmicro.f \
		champcirculaire.f \
		champlineaire.f \
		champlineairemicro.f \
		champlineairemulti.f \
		champlineairemultimicro.f \
		champmultifft.f \
		champtotalcirc.f \
		champtotalgauss.f \
		champtotallineaire.f \
		champtotallineairemulti.f \
		coeffmie.f \
		comparaison.f \
		comparaisonxyz.f \
		computegcfft2dsurf.f \
		d1mach.f \
		d9lgmc.f \
		dasyjy.f \
		dbesj.f \
		dcsevl.f \
		derivative.f \
		dgamlm.f \
		dgamma.f \
		diffractefft2dsurf2.f \
		diffractefft2dsurf2lens.f \
		diffractefft2dtoepos.f \
		djairy.f \
		dlngam.f \
		dznrm2.f \
		espace_libreIcc.f \
		fdump.f \
		fonctiongreensurfcomp.f \
		fonctiongreensurfcompinterp.f \
		fonctiongreensurffft.f \
		fonctiongreensurfinterpfft.f \
		gaussianiso.f \
		gaussiansurfcirc.f \
		gaussiansurf.f \
		gaussiansurfmicro.f \
		gausskronrodpattersonmulti.f \
		gpbicg2.f \
		gpbicgar2.f \
		gpbicgar.f \
                gpbicg.f \
                gpbicgsafe.f \
          	gpbicgplus.f \
	  	gpbicor.f \
                cors.f \
		i1mach.f \
		incidentarbitrary.f \
		initds.f \
		interpdielec5.f \
		inversemat33.f \
                inversemat33r.f \
                inverserigsurf.f \
                inverserigsurfopt.f \
		irradiancesurf.f \
		j4save.f \
		local-macro-surf.f \
		mainsurf.f \
		module_mie.f \
		numerocouche.f \
		objectarbitrarysurf.f \
		objectcubesurf.f \
		objectcylindresurf.f \
		objectellipsesurf.f \
		objectnspheressurf.f \
		objectparasurf.f \
                objectparanxnynzsurf.f \
                objectsphereconcentricsurf.f \
                objectspheresurf.f \
                objectsphereinhomosurf.f \
                objectparainhomosurf.f \
                objectparanxnynzinhomosurf.f \
		passagefourierimage.f \
		passagefourierimagegross.f \
		pimzbicgstab.f \
		polaepstenscomp.f \
		polarisabilitecomp.f \
                produitfftmatvectsurm.f \
                produitfftmatvectsurmopt.f \
                produitfftmatvectsurmtrans.f \
                produitfftmatvectsurmopttrans.f \
		produitfftmatvectsurplusboite.f \
		produitfftmatvectsurplus.f \
		qmrbicgstab.f \
		qmrpim.f \
		ratint.f \
		relecturesurf.f \
		sommearg.f \
		specklesurf.f \
		splinecomp.f \
		tenseurmulticouchec.f \
		tenseurmulticoucheloin.f \
		tenseurmulticoucheloinfft.f \
		tfqmr.f \
		xerbla.f \
		xercnt.f \
		xerhlt.f \
		xermsg.f \
		xerprn.f \
		xersve.f \
		xgetua.f \
		zcg.f \
                zgedid.f \
                champlineairekxky.f \
                champlineairemicrokxky.f \
                microssurfbf.f \
                microssurfdf.f \
                passagefourierimage2.f \
                passagefourierimagegross2.f \
                diffractefft2dtoeposfour.f \
                deltakroutine.f \
		polamodifie.f
                
INCLUDEPATH 	+= .

LIBS 		+= 	-lgfortran

QMAKE_DISTCLEAN += lib/*
