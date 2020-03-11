uniform float randseed1;
uniform float randseed2;

void main() {
    float x = gl_FragCoord.x;
	float y = gl_FragCoord.y;
	
	if (randseed1 > randseed2) {
		x = randseed2 * sin(x);
	} else {
		y = randseed1 * sin(y);
	}
	
	if (sin(x + y) > sin(x - y * x)) {
		x = x + y * randseed1;
	} else {
		x = x - y * randseed2;
	}
	
	float c1 = sin(x * randseed1 + y);
	float c2 = sin(x * y * randseed2 + x);
	
	float scale;
	if (c1 > c2) {
		scale = 1;
	} else {
		scale = 0;
	}
	
	float c3 = cos(c1 + c2 * c1 * x);
	float c4 = cos(c2 - c1 * c1);
	
	if (c3 > c4) {
		scale = 1 - scale;
	}
	
	gl_FragColor.x = scale;
	gl_FragColor.y = scale;
	gl_FragColor.z = scale;
	gl_FragColor.w = 1;
}
