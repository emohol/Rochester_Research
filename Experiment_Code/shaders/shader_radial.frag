#version 120

// Do not touch the code bellow
#define M_PI 3.1415926535897932384626433832795
uniform float screenOffsetX;	// Transfrom coordinate system into Eyeris
uniform float screenOffsetY;	// Transfrom coordinate system into Eyeris
uniform float cX;	// Center of object
uniform float cY;	// Center of object

// Add your code here
uniform float SpatialFreq;
uniform float eccentricity;
uniform float width;
uniform float gaussSD; // std of gaussian envelope of annulus
uniform float gaussSD2; // std of gaussian envelope for central grating
uniform float pixelAngle;
uniform float amplitude; // contrast scalar
uniform float transparency; 
uniform float phase;

void main(void)
{
	// Get coordinate 
	float x = (gl_FragCoord.x - screenOffsetX) * pixelAngle / 60;
	float y = (gl_FragCoord.y - screenOffsetY) * pixelAngle / 60; // in degrees

	float r = sqrt(x * x + y * y);
	
	float scale;

	if (eccentricity > 0) {
		if (r >= eccentricity - width/2 & r <= eccentricity+width/2) {
			scale = (amplitude * sin(2 * M_PI * SpatialFreq * r + phase) + 1) / 2;
		} else if (r > eccentricity + width/2) {
			scale = (amplitude * sin(2 * M_PI * SpatialFreq * r  + phase) * 
					exp(-pow(r - (eccentricity + width/2), 2) / 2 / pow(gaussSD, 2)) + 1) / 2;
		} else if (r < eccentricity - width/2) {
			scale = (amplitude * sin(2 * M_PI * SpatialFreq * r  + phase) * 
					exp(-pow(r - (eccentricity - width/2), 2) / 2 / pow(gaussSD, 2)) + 1) / 2;
		}
	} else {
		scale = (amplitude * sin(2 * M_PI * SpatialFreq * r  + phase) * 
			exp(-pow(r, 2) / 2/ pow(gaussSD2, 2)) + 1) / 2;
	}
	
	gl_FragColor.x = scale;
	gl_FragColor.y = scale;
	gl_FragColor.z = scale;
	gl_FragColor.w = transparency;
}
