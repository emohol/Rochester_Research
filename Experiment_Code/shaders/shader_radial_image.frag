#version 130

// Do not touch the code bellow
#define M_PI 3.1415926535897932384626433832795
uniform sampler2D texture_0; // Texture index 
uniform float screenOffsetX;	// Transfrom coordinate system into Eyeris
uniform float screenOffsetY;	// Transfrom coordinate system into Eyeris
uniform float cX;	// Center of object
uniform float cY;	// Center of object

// Add your code here
uniform float rcx;
uniform float rcy;
uniform float SpatialFreq;
uniform float eccentricity;
uniform float width;
uniform float gaussSD; // std of gaussian envelope of annulus
uniform float gaussSD2; // std of gaussian envelope for central grating
uniform float pixelAngle;
uniform float amplitude; // contrast scalar
uniform float transparency; 
uniform float phase;
uniform float imageScale; // contrast of input image

void main(void)
{
	// Get the coordinate 
	float x = gl_FragCoord.x - (screenOffsetX + cX);
	float y = gl_FragCoord.y - (screenOffsetY + cY);
	
	float r = sqrt((x - rcx)*(x-rcx) + (y-rcy) * (y-rcy)) * pixelAngle / 60;
	
	float scale;

	// Raised cosine envelope (JI 11-26-2018). In this case let gaussSD be half
	// the period of the cosine
	if (eccentricity == 0) {
		//gaussSD = gaussSD2; // JI: we probably don't need this option
	}
	
	if ((r >= eccentricity - width/2) && (r <= eccentricity+width/2)) {
		scale = amplitude * sin(2 * M_PI * SpatialFreq * r + phase);
	} else if (r > eccentricity + width/2 && r <= eccentricity + width/2 + gaussSD) {
		scale = amplitude * sin(2 * M_PI * SpatialFreq * r  + phase) * 
					(cos(M_PI / gaussSD * (r - eccentricity - width / 2)) + 1) / 2;
	} else if (r < eccentricity - width/2 && r >= eccentricity - width/2 - gaussSD) {
		scale = amplitude * sin(2 * M_PI * SpatialFreq * r  + phase) * 
					(cos(M_PI / gaussSD * (r - eccentricity + width / 2)) + 1) / 2;
	} else {
		scale = 0;
	}

	
	// Gaussian envelope
	/*s
	if (eccentricity > 0) {
		if (r >= eccentricity - width/2 && r <= eccentricity+width/2) {
			scale = amplitude * sin(2 * M_PI * SpatialFreq * r + phase);
		} else if (r > eccentricity + width/2) {
			scale = amplitude * sin(2 * M_PI * SpatialFreq * r  + phase) * 
					exp(-pow(r - (eccentricity + width/2), 2) / 2 / pow(gaussSD, 2));
		} else if (r < eccentricity - width/2) {
			scale = amplitude * sin(2 * M_PI * SpatialFreq * r  + phase) * 
					exp(-pow(r - (eccentricity - width/2), 2) / 2 / pow(gaussSD, 2));
		}
	} else {
		scale = amplitude * sin(2 * M_PI * SpatialFreq * r  + phase) * 
			exp(-pow(r, 2) / 2/ pow(gaussSD2, 2));
	}
	*/
	
	// Get color from the texture
	vec4 Color = texture2D(texture_0, vec2(gl_TexCoord[0]));
	
	// Modify color of object by the position of pixel
	gl_FragColor.x = (imageScale * (Color.x - .5) + .5) + scale;
	gl_FragColor.y = (imageScale * (Color.y - .5) + .5) + scale;
	gl_FragColor.z = (imageScale * (Color.z - .5) + .5)+ scale;
	gl_FragColor.a = Color.a;
}