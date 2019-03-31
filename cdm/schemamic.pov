// reseau.pov
// a lancer avec la commande  povray +Q11 +A0.1 +R8 +P reseau.pov pour un mailleur rendue

#include "colors.inc"
#include "textures.inc"
#include "shapes.inc"
#include "metals.inc"
#include "glass.inc"
#include "woods.inc"   
#fopen totof "toto" write



background{White}

global_settings {
//    assumed_gamma 1
    max_trace_level 5
//    photons {
//        spacing 0.00025
//        count 150000
//        max_trace_level 9
//        media 100, 2
//        media 500, 3
//ambient_light <1,1,1>
    }

   
/* global_settings { radiosity {
	distance_maximum 0.2
	nearest_count 10
	count 1000
	gray_threshold 0.3
}} */
   
   
#declare CamPos = < -25, 65, -110>;

camera {
    location CamPos
    look_at < 25, -1, 0>
    angle 35
}

light_source {CamPos, color Gray25
    photons {refraction off reflection off}
    media_interaction off
}


light_source  { < 200 , 300 , -200 > color <1,1,1> }

light_source  { <-600, 500, 400> color <0.3,0.3,0.3>}

light_source  { < -500 , 300 , -200 > color <1,1,1> }

//light_source  { < 0 , -50 , 0 > color <1,1,1> }


//--------------------------  lens --------
 #declare zlens1 = 20; //position lens			      
 #declare R   = 60.0;   //sphere radius
 #declare Over= 1;   //sphere overlap
 intersection{
  sphere{<0,0,0>,R  translate <zlens1-R+Over,0,0>}
  sphere{<0,0,0>,R  translate <zlens1+R-Over,0,0 >}
  texture{ Dark_Green_Glass }
  interior {ior 1.5} 
 finish {diffuse 5
                  brilliance 0.0
                  specular 0.3
                  roughness 0.2
                  //reflection 0.5 
                  refraction 0. } 
  translate < 0,0,0>}


//--------------------------  lens --------
 #declare zlens2 = 40; //position lens			      
 #declare R   = 80.0;   //sphere radius
 #declare Over= 1;   //sphere overlap
 intersection{
  sphere{<0,0,0>,R  translate <zlens2-R+Over,0,0>}
  sphere{<0,0,0>,R  translate <zlens2+R-Over,0,0 >}
  texture{ Dark_Green_Glass }
  interior {ior 1.5} 
 finish {diffuse 5.
                  brilliance 10.0
                  specular 0.3
                  roughness 0.2
                  //reflection 0.5 
                  refraction 0. } 
  translate < 0,0,0>}
 //---------------------------------- end -- 


#declare ObjCDM = union {
#declare a=1;
#declare k = -3 ;
#while ( k < 4 ) // la condition
        #declare j = -3 ; 
	#while ( j < 4 ) 
        	#declare i = -3 ;
        	#while ( i < 4 ) // la condition
			#if (i*i+j*j+k*k < 10)
				sphere { <i, j, k> 0.8       // radius of sphere
				texture {pigment {color < 0.4 , 0.4 , 1 >  filter .9} } finish {Glossy}
				}
			//	box { <i-0.5, j-0.5, k-0.5> <i+0.5, j+0.5, k+0.5>
			//	texture {pigment {color < 0.4 , 0.4 , 1 >  filter .9}} finish {Glossy}
			//	}
			#end
	        #declare i = i + 1 ; 
        	#end                 
	#declare j = j + 1 ; 
	#end                 
#declare k = k + 1 ; 
#end                 
}

object {ObjCDM
scale <1.5,1.5,1.5> translate  <7,0,0>
}



#declare len=1;
#declare epsilon=0.2;
#declare cote=12;

#declare CCD =
box {  < 60 , -cote  , -cote >    < 60+epsilon , cote , cote >
//pigment {color < .6 , .3 , 0.1 >}
 texture { T_Glass4 } interior {ior 1.5} 
 finish {diffuse 0.
                  brilliance 0.0
                  specular 0.3
                  roughness 0.2
                 reflection 0.1
                  refraction 0.1 } 
}     


object {CCD
scale <1,1,1>
}

           

//DEFINITION DES AXES

//AXE X
#declare length=15;
#declare rayon=0.1;
#declare cylindre = cylinder{<0,0,0>,<0,0,length>, rayon rotate<0,0,0>
translate<0,0,0> texture{pigment{color Black}
                  finish {ambient 0.15 diffuse 0.75 reflection 0.1 phong 1}}}
#declare bout =cone{ <0,0,length> , rayon+0.2 , <0,0,length+1> , 0
   translate<0,0,0> texture{pigment{color Black}
   finish {ambient 0.15 diffuse 0.75 reflection 0.1 phong 1}}}        
#declare fleche = merge{ object{cylindre} object {bout}}
object{fleche rotate<0,0,0> translate<0,0,0> }


//AXE Z
#declare length=60;
#declare lengthm=10;
#declare rayon=0.1;
#declare cylindre = cylinder{<0,0,-lengthm>,<0,0,length>, rayon rotate<0,0,0>
translate<0,0,0> texture{pigment{color Black}
                  finish {ambient 0.15 diffuse 0.75 reflection 0.1 phong 1}}}
#declare bout =cone{ <0,0,length> , rayon+0.2 , <0,0,length+1> , 0
   translate<0,0,0> texture{pigment{color Black}
   finish {ambient 0.15 diffuse 0.75 reflection 0.1 phong 1}}}        
#declare fleche = merge{ object{cylindre} object {bout}}
object{fleche rotate<0,90,0> translate<0,0,0> }

//AXE Y
#declare length=13;
#declare rayon=0.1;
#declare cylindre = cylinder{<0,0,0>,<0,0,length>, rayon rotate<0,0,0>
translate<0,0,0> texture{pigment{color Black}
                  finish {ambient 0.15 diffuse 0.75 reflection 0.1 phong 1}}}
#declare bout =cone{ <0,0,length> , rayon+0.2 , <0,0,length+1> , 0
   translate<0,0,0> texture{pigment{color Black}
   finish {ambient 0.15 diffuse 0.75 reflection 0.1 phong 1}}}        
#declare fleche = merge{ object{cylindre} object {bout}}
object{fleche rotate<-90,0,0> translate<0,0,0> }


//DEFINITION ILLUMINATION




#declare length=15;
#declare rayon=0.3;
#declare cylindre = cylinder{<0,0,-length>,<0,0,-1>, rayon rotate<0,0,0>
translate<0,0,0> texture{pigment{color Green}
                  finish {ambient 0.15 diffuse 0.75 reflection 0.1 phong 1}}}
#declare bout =cone{ <0,0,-1> , rayon+0.2 , <0,0,0> , 0
   translate<0,0,0> texture{pigment{color Green}
   finish {ambient 0.15 diffuse 0.75 reflection 0.1 phong 1}}}        
#declare fleche1 = merge{ object{cylindre} object {bout}}
object{fleche1 rotate<0,45,0> translate<0,0,0> }

//TE

#declare length=10;
#declare rayon=0.3;
#declare cylindre = cylinder{<0,0,-length>,<0,0,-1>, rayon rotate<0,0,0>
translate<0,0,0> texture{pigment{color Red}
                  finish {ambient 0.15 diffuse 0.75 reflection 0.1 phong 1}}}
#declare bout =cone{ <0,0,-1> , rayon+0.2 , <0,0,0> , 0
   translate<0,0,0> texture{pigment{color Red}
   finish {ambient 0.15 diffuse 0.75 reflection 0.1 phong 1}}}        
#declare fleche1 = merge{ object{cylindre} object {bout}}
object{fleche1 rotate<-90,0,0> translate<-10.6,length,-10.6> }


text { ttf "timrom.ttf" "TE//(x,y)" 0.05, 0 pigment { Red }
scale<2,2,2> rotate<0,0,0> translate<-13.6,length+1, -10.6> }


//TM

#declare length=10;
#declare rayon=0.3;
#declare cylindre = cylinder{<0,0,length>,<0,0,0>, rayon rotate<0,0,0>
translate<0,0,0> texture{pigment{color Red}
                  finish {ambient 0.15 diffuse 0.75 reflection 0.1 phong 1}}}
#declare bout =cone{ <0,0,length> , rayon+0.2 , <0,0,length+1> , 0
   translate<0,0,0> texture{pigment{color Red}
   finish {ambient 0.15 diffuse 0.75 reflection 0.1 phong 1}}}        
#declare fleche1 = merge{ object{cylindre} object {bout}}
object{fleche1 rotate<180,-45,0> translate<-10.6,0,-10.6> }


text { ttf "timrom.ttf" "TM" 0.05, 0 pigment { Red }
scale<2,2,2> rotate<0,0,0> translate<-3,-1, -20.6> }


text { ttf "timrom.ttf" "z" 0.05, 0 pigment { Black }
scale<2,2,2> translate<9.8,0, -1> }

text { ttf "timrom.ttf" "x" 0.05, 0 pigment { Black }
scale<2,2,2> translate<-2.0,0.5,11> }

text { ttf "timrom.ttf" "y" 0.05, 0 pigment { Black }
scale<2,2,2> translate<1,12.0,0> }

//text { ttf "symbol.ttf" "Q" 0.05, 0 pigment { Green }
//scale<3,3,3> translate<-12,0,0> rotate<0,-25,0>}

//text { ttf "symbol.ttf" "j" 0.05, 0 pigment { Green }
//scale<3,3,3> rotate<0,45,0> translate<0,-10,0>  }

//Kinc

text { ttf "timrom.ttf" "k" 0.05, 0 pigment { Green }
scale<3,3,3> translate<-20.5,3.0, 0>rotate<0,-45,0>  }

text { ttf "timrom.ttf" "inc" 0.05, 0 pigment { Green }
scale<2,2,2> translate<-18.5,2.0, 0>rotate<0,-45,0>  }



// tore

#declare  rayext= 7;
#declare  rayint= 0.2;
#declare tora = torus{ rayext rayint pigment { Green }}
#declare boite1 =box { <0,-2*rayext,-2*rayext > <2*rayext,2*rayext,2*rayext > pigment { Green }}
#declare boite2 =box { <0,-2*rayext,-2*rayext > <2*rayext,2*rayext,2*rayext >  rotate<0,135,0> pigment { Green } }
#declare demitore = difference{ object{tora} object{boite1}  }
#declare quarttore = difference{ object{demitore} object{boite2}  }

//object{quarttore rotate<0,-90,0> translate<0,0,0>}


text { ttf "timrom.ttf" "Objective lens" 0.05, 0 pigment { Black }
scale<3,3,3>  rotate<0,90,0>  translate<20,15.0, 0>  }

text { ttf "timrom.ttf" "Tube lens" 0.05, 0 pigment { Black }
scale<3,3,3> rotate<0,90,0>  translate<40,15.0, 0>  }

text { ttf "timrom.ttf" "CCD camera:" 0.05, 0 pigment { Black }
scale<3,3,3> rotate<0,90,0>  translate<60,18.0, 0>  }
text { ttf "timrom.ttf" "Image plane" 0.05, 0 pigment { Black }
scale<3,3,3> rotate<0,90,0>  translate<60,15.0, 0>  }

//text { ttf "timrom.ttf" "N=7" 0.05, 0 pigment { Black }
//scale<3,3,3> rotate<0,0,0>  translate<5,-10,0>  }
text { ttf "timrom.ttf" "Discretization" 0.05, 0 pigment { Black }
scale<3,3,3> rotate<0,0,0>  translate<3,-13.0,0>  }



text { ttf "timrom.ttf" "IF-DDA" 0.5, 0 pigment { Red }
scale<10,10,10> rotate<0,0,0>  translate<5,-30,0>  }
