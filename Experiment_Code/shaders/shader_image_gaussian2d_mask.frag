// Do not touch the code bellow
uniform sampler2D texture_0; // Texture index 
uniform float screenOffsetX;	// Transfrom coordinate system into Eyeris
uniform float screenOffsetY;	// Transfrom coordinate system into Eyeris
uniform float cX;	// Center of object
uniform float cY;	// Center of object

// Add your code here
uniform float sigma;
uniform float transparency;

void main(void)
{
	// Get the coordinate 
	float x = gl_FragCoord.x - (screenOffsetX + cX);
	float y = gl_FragCoord.y - (screenOffsetY + cY);
	float r2 = x * x + y * y;
	float sigma2 = 2 * sigma * sigma;
	
	float scale = 0;
	if (sigma2 > 0) {
		scale = exp( - r2 / sigma2);
	}
	
	// Get color from the texture
	vec4 Color = texture2D(texture_0, vec2(gl_TexCoord[0]));
	
	// Modify color of object by the position of pixel
	gl_FragColor.w = transparency;
	gl_FragColor.x = Color.x * scale;
	gl_FragColor.y = Color.y * scale;
	gl_FragColor.z = Color.z * scale;
}