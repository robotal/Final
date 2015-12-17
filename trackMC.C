#include <TMath.h>
#include <TRandom.h>
#include <algorithm>
//#include <algorithm.h>
//include header file algorithm

// parameters that vould be modified
//
Int_t n = 10000; //number of trials
Float_t momentum = 1; //GeV
Float_t field = 3; //Tesla
Float_t maxRadius = 1.5; //m
static const Int_t layers = 12;
Float_t charge = 1; //e
Bool_t plotReciprocals = false; //true if need to plot reciprocals 

//things we define below should not be modified

// a function to get line of best fit
float * lineOfBestFit(float x[], float y[], float err[], float size){
	float s1_err2 = 0; //sum of uncertainty^2
	float sx_err2 = 0; //sum of x/uncertainty^2
	float sx2_err2 = 0; //sum of x^2/uncertainty^2
	float sy_err2 = 0; //sum of y/uncertainty^2
	float sxy_err2 = 0; //sum of xy/uncertainty^2

	for(int i = 0; i < size; i++){

		s1_err2 += 1/pow(err[i],2);
		sx_err2 += x[i]/pow(err[i],2);
		sx2_err2 += pow(x[i],2)/pow(err[i],2);
		sy_err2 += y[i]/pow(err[i],2);
		sxy_err2 += x[i]*y[i]/pow(err[i],2);

	}

	float delta = s1_err2*sx2_err2-sx_err2*sx_err2; //delta factor
	float b = 1/delta*(sx2_err2*sy_err2-sx_err2*sxy_err2); //y intercept
	float m = 1/delta*(s1_err2*sxy_err2-sx_err2*sy_err2); //slope

	static float line[2];
	line[0]=m; //slope
	line[1]=b; //y intercept


	return line;
}

// function that returns Absolute value of a number
float absolute(float n){

	if(n>=0){

		return n;
	}
	return n*-1.0;
}

/*
 * Here is what is going down in each iteration of the simulation:
 * 1) Calculate the actual radius of the particle's path (pT = .3qrB)
 * 2) Generate a and b, the coordinates of the center of the circle, pseudo-randomly s.t. a^2+b^2 = r^2, where r is the radius
 * of the particle's path
 * 3) Now, at each layer of the tracker:
 * a) Calculate the intersection of the circle representing the layer and the circle of the particle's actual path
 * b) Convert those x-y coordinates to r-phi coordinates
 * c) Smear phi up to 10 micrometers in the direction orthogonal to the radial segment to the origin
 * d) Convert back to x-y (these are the measured x-y coordinates)
 * e) Add the u-v coordinates to an array (linearization of the circle)
 * 4) Fit the u-v coordinates to generate the equation of the measured circle
 * 5) Calculate the measured radius of the particle's path based on this equation (a^2+b^2 = r^2)
 * 6) Compute and plot the measured momentum of the particle (pT = .3qrB)

*/

void trackMC(){
	//================================================
	//Set up the display
	//================================================

	TCanvas *c1 = new TCanvas("c1","Tracker Simulation",200,10,700,500); 
	c1->SetFillColor(42);

	//calculate 1/2 the range of the histogram
	float range = .000012*pow(momentum,2);



	TH1F *histo1;
	//define two histograms 
	//if plotReciprocals is true then define histo1 to be as Reciprocals of Measured PTs
	if(plotReciprocals){

		histo1 = new TH1F("histo1","Reciprocals of Measured PTs",100,1/(momentum+range),1/(momentum-range));
	}else{

		histo1 = new TH1F("histo1","Measured PTs",100,momentum-range,momentum+range);

	}

	histo1->SetMarkerStyle(21);
	//================================================
	// random number generation
	//================================================

	// define TRandom ran 
	//gRandom get function SetSeed()

	TRandom ran;
	gRandom->SetSeed();

	//================================================
	// Compute the actual radius of the particle's path (pT = .3qrB)
	//================================================

	//define radius as float radiusActual = pT /.3qB

	Float_t radiusActual = momentum / (.3 * charge * field);	

	//================================================
	// Compute the actual radius of the particle's path (pT = .3qrB)
	//================================================

	for(int i = 0; i < n; i++){
		//randomly generate a and b such that a^2 + b^2 = r^2, where r is the radius of curvature of the particle
		//for simplicity, a and b lie in the first quadrant
		float a = ran.Uniform(0,radiusActual);
		float b = sqrt(pow(radiusActual,2)-pow(a,2));

		//u-v coordinates to be later used for tracking
		float u[layers*2];
		float v[layers*2];
		float err[layers*2];


		//take measurements at each layer
		//start a for loop to loop over number of layers
		for(int j = 1; j <= layers; j++){
			float r = maxRadius/layers*j; //current radius


			
			//A, B, and C for quadratic equation Ax^2+Bx+C = 0, used to determine intersection of the circle
			//representing the current radius (x^2+y^2 = r^2) and the circle (x-a)^2+(y-b)^2=a^2+b^2;
			float qA = 1+pow(a,2)/pow(b,2);
			float qB = -a*pow(r,2)/pow(b,2);
			float qC = pow(r,2)*(pow(r,2)/(4*pow(b,2))-1);
			float x1 = (-qB+sqrt(pow(qB,2)-4*qA*qC))/(2*qA); //The x coordinate of the intersection (the greater one)
			float x2 = (-qB-sqrt(pow(qB,2)-4*qA*qC))/(2*qA); //The x coordinate of the intersection (the lesser one)

			float y1 = -a/b*x1+pow(r,2)/(2*b); //The y-coordinate coresponding to the greater x coordinate (x1)
			float y2 = -a/b*x2+pow(r,2)/(2*b); //The y-coordinate coresponding to the lesser x coordinate (x2)

			//Compute r-phi coordinates
			float phi1, phi2;
			if(x1 > 0){
				phi1 = atan(y1/x1);
			}
			if(x1 < 0){
				phi1 = atan(y1/x1)+TMath::Pi(); 
			}
			if(x2 > 0){
				phi2 = atan(y2/x2);
			}
			if(x2 < 0){
				phi2 = atan(y2/x2)+TMath::Pi(); 
			}


			//smear the coordinates of x1-y1
			float epsilon = 10*pow(10,-6); //error in measurement perpendicular to the radial segment
			float Ephimax = atan(epsilon/r); //maximum error in phi

			//"smear" the phi value pseudo-randomly according to a gaussian distribution
			//define float efactor equal to random Gaus(0,1)

			Float_t efactor = ran.Gaus(0,1);

			if(efactor > 1){
				efactor = 1;
			}
			if(efactor < -1){
				efactor = -1;
			}
			phi1 += Ephimax*efactor;

			//convert from r-phi back to x-y
			float X1m = r*cos(phi1);
			float Y1m = r*sin(phi1);

			//take the maximum uncertainty in y between two calculations, one of which assumes
			//that theta is its maximum positive value and the other of which assumes that theta
			//is its maximum negative value

			float erry1 = TMath::max(absolute(epsilon*cos(Ephimax)/cos(phi1)), absolute(epsilon/cos(phi1-Ephimax)));

			//smear the coordinates of x2-y2
			// //set efactor again to random Gaus(0,1)
			
			efactor = ran.Gaus(0,1);
			
			if(efactor > 1){
				efactor = 1;
			}
			if(efactor < -1){
				efactor = -1;
			}
			phi2 += Ephimax*efactor;
			float X2m = r*cos(phi2);
			float Y2m = r*sin(phi2);
			float erry2 = TMath::max(absolute(epsilon*cos(Ephimax)/cos(phi2)), absolute(epsilon/cos(phi2-Ephimax)));

			//add values to u-v and error arrays
			u[2*(j-1)] = X1m/(r*r);
			v[2*(j-1)] = Y1m/(r*r);

			err[2*(j-1)] = erry1/(r*r);

			u[2*(j-1)+1] = X2m/(r*r);
			v[2*(j-1)+1] = Y2m/(r*r);
			err[2*(j-1)+1] = erry2/(r*r);
		}//end for loop over number of layers


		//calculate a and b as measured by the tracker after fitting
		float *line = lineOfBestFit(u,v, err, layers*2);

		float slope = *(line+0);
		float yInt = *(line+1);

		//measured a and b values
		float bm = 1/(2*yInt);
		float am = -slope*bm;


		//calculate the radius measured by the tracker (a^2+b^2 = r^2)
		float rm = sqrt(pow(am,2)+pow(bm,2));

		//calculate the energy measured by the tracker (p = .3qrB)
		float pm = .3*charge*field*rm;

		//plot the momentum or reciprocal of the momentum
		//if plotReciprocals true fill histo1 reciprocal of pm
		//else fill pm
		//
		//

		if(plotReciprocals){

			cout << 1/pm << endl;
			histo1->Fill(1/pm);
		}
		else{

			cout << pm << endl;
			histo1->Fill(pm);
		}



	}//end of loop on n

	histo1->Fit("gaus");
	//Fit the histogram to a gaussian distribution
	histo1->SetLineColor(2);
	histo1->Draw("");
	//update the canvas
	c1->Update();

}
