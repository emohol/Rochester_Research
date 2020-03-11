// Do not touch the code bellow
uniform sampler2D texture_0; // Texture index 
uniform float screenOffsetX;	// Transfrom coordinate system into Eyeris
uniform float screenOffsetY;	// Transfrom coordinate system into Eyeris
uniform float cX;	// Center of object
uniform float cY;	// Center of object

// Add your code here
#define M_PI 3.1415926535897932384626433832795
uniform float angle;
uniform float wavelength;
uniform float transparency;

void main(void)
{
	// Get coordinate 
	float x = gl_FragCoord.x - (screenOffsetX + cX);
	float y = gl_FragCoord.y - (screenOffsetY + cY);
	
	// Calculate grating
	float phase = (x * cos(angle) + y * sin(angle)) / wavelength * 2 * M_PI;
	float scale = (cos(phase) + 1) / 2;
	
	// Get color from the image
	vec4 Color = texture2D(texture_0, vec2(gl_TexCoord[0]));
	
	// Modify color of object by the position of pixel
	gl_FragColor.w = transparency;
	gl_FragColor.x = Color.x * scale;
	gl_FragColor.y = Color.y * scale;
	gl_FragColor.z = Color.z * scale;

}