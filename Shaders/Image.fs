uniform sampler2D iChannel0;
uniform sampler2D iChannel1;
uniform sampler2D iChannel2;
uniform sampler2D iChannel3;
uniform float iTime;
uniform vec4 iMouse;
uniform vec2 iResolution;
out vec4 C;
#define MAX_STEPS 100
#define MAX_DIST 100.
#define SURF_DIST .01


float getDist(vec3 p)
{
    vec4 sphere = vec4(0.0,1.0,6.0,1.0);
    float sphereDist = length(p - sphere.xyz) - sphere.w;
    float planeDist = p.y;
    float d = min(sphereDist,planeDist);
    return d;

}



void main()
{
    vec2 p = (gl_FragCoord.xy-0.5*iResolution)/iResolution.y;


    vec3 col = mix(vec3 (.01,0.1,0.1),vec3(0.8,0.8,1.0),-p.y);
    vec3 Can = vec3(0.,1.0,-1.0);
    vec3 dir = normali


    C = vec4( col, 1.0 );
}

