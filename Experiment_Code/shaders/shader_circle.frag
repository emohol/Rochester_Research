// Do not touch the code bellow
uniform float screenOffsetX;	// Transfrom coordinate system into Eyeris
uniform float screenOffsetY;	// Transfrom coordinate system into Eyeris
uniform float cX;	// Center of object
uniform float cY;	// Center of object

// Add your code here
#define M_PI 3.1415926535897932384626433832795
uniform float wavelength;
uniform float Radius;
uniform float transparency;

void main(void)
{
	// Get coordinate 
	float x = gl_FragCoord.x - (screenOffsetX + cX);
	float y = gl_FragCoord.y - (screenOffsetY + cY);
	
	// Calculate circular grating
	float r = sqrt(dx * dx + dy * dy);
	float scale = (cos(r / wavelength * 2 * M_PI) + 1) / 2;

	// Modify the color of object by the position of pixel
	if (r < Radius) {
		gl_FragColor.w = transparency;
		gl_FragColor.x = gl_Color.x * scale;
		gl_FragColor.y = gl_Color.y * scale;
		gl_FragColor.z = gl_Color.z * scale;
	} else {
		gl_FragColor = vec4(0, 0, 0, 1);
	}
}
