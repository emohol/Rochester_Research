#version 130

// Do not touch this file

void main(void)
{
    gl_TexCoord[0] = gl_MultiTexCoord0;
	gl_FrontColor = gl_Color;
    gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
}
