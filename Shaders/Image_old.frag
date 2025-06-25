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
#define EPSILON 0.01

float sphereDist(vec3 p, vec4 sphere)
{
    return length(p-sphere.xyz)-sphere.w;
}
float planeDist(vec3 p)
{
    return p.y;
}

const vec3 COLOR_BACKGROUND = vec3(0.35, 0.1, 0.0);

float getDist(vec3 p)
{
    vec4 sphere = vec4(0.0,1.0,-5.0,1.0);
    float sphDist = sphereDist(p,sphere);
    float plaDist = planeDist(p);
    float d = min(sphDist,plaDist);
    return d;

}

float RayMarching(vec3 Cam,vec3 dir)
{
    int i=0;
    float d0=0.0;
    vec3 p = Cam+d0*dir;
    float dist = getDist(p);
    while((i<MAX_STEPS)&&(dist>SURF_DIST)&&(d0<MAX_DIST))
    {
        d0+=dist;
        p= Cam+d0*dir;
        dist= getDist(p);
        i++;
    }
    return d0;
}

vec3 estimateNormal(vec3 p)
{
    float d= getDist(p);

    float dx = d-getDist(vec3(p.x-EPSILON,p.y,p.z));
    float dy = d-getDist(vec3(p.x,p.y-EPSILON,p.z));
    float dz = d-getDist(vec3(p.x,p.y,p.z-EPSILON));
    vec3 n =vec3(dx,dy,dz);
    return normalize(n);
}

float getLight(vec3 p)
{
    vec3 lightPos = vec3(0.0,5.0,3.0);
    lightPos.xy+=vec2(cos(iTime),sin(iTime))*2.0;
    vec3 n = estimateNormal(p);
    vec3 lightDir = normalize(lightPos-p);
    float diff = clamp(dot(n,lightDir),0.0,1.0);
    float d = RayMarching(p+n*SURF_DIST, lightDir);
        if(d<length(lightPos-p)) diff *= .1;
    return diff;
}

void main()
{
    vec2 p = (gl_FragCoord.xy-0.5*iResolution)/iResolution.y;


    vec3 col = mix(vec3 (.01,0.1,0.1),COLOR_BACKGROUND,-p.y*0.2);
    vec3 Can = vec3(0.0,1.0,2.0);
    float fieldOfView = 90.0;
    float z = (iResolution.x/iResolution.y) / tan(radians(fieldOfView) / 2.0);
    vec3 dir = normalize(vec3(p.x,p.y,-z));

    float d= RayMarching(Can,dir);
    vec3 ps = Can+d*dir;
    float diff = getLight(ps);
    //d/=10;
    col=vec3(diff);
    col = mix(COLOR_BACKGROUND,col, smoothstep(0.1,0.9,diff));//smoothstep(col,COLOR_BACKGROUND,vec3(d));
    C = vec4( col, 1.0 );
}

