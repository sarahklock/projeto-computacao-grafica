in vec3 aColor;
in vec4 aPosition;
out vec4 C;
uniform sampler2D iChannel0;
uniform sampler2D iChannel1;
uniform sampler2D iChannel2;
uniform sampler2D iChannel3;
uniform vec2 iResolution;
uniform vec4 iMouse;
uniform float iTime;
uniform int iFrame;



#define MAX_STEPS 100
#define MAX_DIST 200.
#define SURF_DIST .01
#define EPSILON .01
#define PI 3.14159265359
float dot2( in vec2 v ) { return dot(v,v); }
float dot2( in vec3 v ) { return dot(v,v); }
float ndot( in vec2 a, in vec2 b ) { return a.x*b.x - a.y*b.y; }

vec3 Sky( vec3 ray )
{
        return mix( vec3(.8), vec3(0), exp2(-(1.0/max(ray.y,.01))*vec3(.4,.6,1.0)) );
}

// noise
float noise(vec2 pos)
{
        return fract( sin( dot(pos*0.001 ,vec2(24.12357, 36.789) ) ) * 12345.123);
}


// blur noise
float smooth_noise(vec2 pos)
{
        return   ( noise(pos + vec2(1,1)) + noise(pos + vec2(1,1)) + noise(pos + vec2(1,1)) + noise(pos + vec2(1,1)) ) / 16.0
                   + ( noise(pos + vec2(1,0)) + noise(pos + vec2(-1,0)) + noise(pos + vec2(0,1)) + noise(pos + vec2(0,-1)) ) / 8.0
           + noise(pos) / 4.0;
}


// linear interpolation
float interpolate_noise(vec2 pos)
{
        float	a, b, c, d;

        a = smooth_noise(floor(pos));
        b = smooth_noise(vec2(floor(pos.x+1.0), floor(pos.y)));
        c = smooth_noise(vec2(floor(pos.x), floor(pos.y+1.0)));
        d = smooth_noise(vec2(floor(pos.x+1.0), floor(pos.y+1.0)));

        a = mix(a, b, fract(pos.x));
        b = mix(c, d, fract(pos.x));
        a = mix(a, b, fract(pos.y));

        return a;
}



float perlin_noise(vec2 pos)
{
        float	n;

        n = interpolate_noise(pos*0.0625)*0.5;
        n += interpolate_noise(pos*0.125)*0.25;
        n += interpolate_noise(pos*0.025)*0.225;
        n += interpolate_noise(pos*0.05)*0.0625;
        n += interpolate_noise(pos)*0.03125;
        return n;
}



const vec3 COLOR_BACKGROUND = vec3(0.25, 0.1, 0.15);

struct Surface
{
    float sd;
    vec3 color;
    float Ka;
    float Kd;
    float Ks;
    int id;
};

mat2 rotate2d(float theta)
{
    float co = cos(theta);
    float s=sin(theta);
    return mat2(co,-s,s,co);
}

mat4 trans (vec3 t)
{
    mat4 mat = mat4 (vec4 (1., .0, .0, .0),
                     vec4 (.0, 1., .0, .0),
                     vec4 (.0, .0, 1., .0),
                     vec4 (t.x, t.y, t.z, 1.));
    return mat;
}

vec3 opTransf (vec3 p, mat4 m)
{
    return vec4 (m * vec4 (p, 1.)).xyz;
}

mat4 rotX (in float angle)
{
    float rad = radians (angle);
    float c = cos (rad);
    float s = sin (rad);

    mat4 mat = mat4 (vec4 (1.0, 0.0, 0.0, 0.0),
                     vec4 (0.0,   c,   s, 0.0),
                     vec4 (0.0,  -s,   c, 0.0),
                     vec4 (0.0, 0.0, 0.0, 1.0));

    return mat;
}

mat4 rotY (in float angle)
{
    float rad = radians (angle);
    float c = cos (rad);
    float s = sin (rad);

    mat4 mat = mat4 (vec4 (  c, 0.0,  -s, 0.0),
                     vec4 (0.0, 1.0, 0.0, 0.0),
                     vec4 (  s, 0.0,   c, 0.0),
                     vec4 (0.0, 0.0, 0.0, 1.0));

    return mat;
}

mat4 rotZ (in float angle)
{
    float rad = radians (angle);
    float c = cos (rad);
    float s = sin (rad);

    mat4 mat = mat4 (vec4 (  c,   s, 0.0, 0.0),
                     vec4 ( -s,   c, 0.0, 0.0),
                     vec4 (0.0, 0.0, 1.0, 0.0),
                     vec4 (0.0, 0.0, 0.0, 1.0));

    return mat;
}


float sphereDist(vec3 p, float r)
{
    vec3 center = p;
    return (length(center)-r*r);
}
float heartDist(vec3 p, vec3 c, float r)
{

    vec3 center = p-c;

    center.x=abs(center.x);
        center.y*=0.9;
        center.z*=(1.1-center.y/15.0);
        center.y+=-0.3-center.x*sqrt(((2.0-center.x)/2.0));
        r+=0.12*pow(0.5+0.5*sin(2.0*PI*iTime+center.y/0.75),4.0);
    return (length(center)-r*r);
}

float planeDist(vec3 p, vec3 Normal,float D)
{
    return p.x*Normal.x + p.y*Normal.y + p.z*Normal.z -D - .1 * sin (4. * p.x) * cos (4. * p.z);
}


float boxDist( vec3 p, vec3 b )
{
    vec3 d = abs(p) - b;
    return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0));
}
float boxRoundDist(vec3 p, vec3 s, float r) {
    vec3 d = abs(p) - (s - r);
    return length(max(d, 0.0)) - r
        + min(max(d.x, max(d.y, d.z)), 0.0);
}


float boxFrameDist( vec3 p, vec3 b, float e )
{
       p = abs(p  )-b;
  vec3 q = abs(p+e)-e;

  return min(min(
      length(max(vec3(p.x,q.y,q.z),0.0))+min(max(p.x,max(q.y,q.z)),0.0),
      length(max(vec3(q.x,p.y,q.z),0.0))+min(max(q.x,max(p.y,q.z)),0.0)),
      length(max(vec3(q.x,q.y,p.z),0.0))+min(max(q.x,max(q.y,p.z)),0.0));
}
float ellipsoidDist( in vec3 p, in vec3 r ) // approximated
{
    float k0 = length(p/r);
    float k1 = length(p/(r*r));
    return k0*(k0-1.0)/k1;
}

float torusDist( vec3 p, vec2 t )
{
    return length( vec2(length(p.xz)-t.x,p.y) )-t.y;
}

float cappedTorusDist(in vec3 p, in vec2 sc, in float ra, in float rb)
{
    p.x = abs(p.x);
    float k = (sc.y*p.x>sc.x*p.y) ? dot(p.xy,sc) : length(p.xy);
    return sqrt( dot(p,p) + ra*ra - 2.0*ra*k ) - rb;
}

float hexPrismDist( vec3 p, vec2 h )
{
    vec3 q = abs(p);

    const vec3 k = vec3(-0.8660254, 0.5, 0.57735);
    p = abs(p);
    p.xy -= 2.0*min(dot(k.xy, p.xy), 0.0)*k.xy;
    vec2 d = vec2(
       length(p.xy - vec2(clamp(p.x, -k.z*h.x, k.z*h.x), h.x))*sign(p.y - h.x),
       p.z-h.y );
    return min(max(d.x,d.y),0.0) + length(max(d,0.0));
}

float octogonPrismDist( in vec3 p, in float r, float h )
{
  const vec3 k = vec3(-0.9238795325,   // sqrt(2+sqrt(2))/2
                       0.3826834323,   // sqrt(2-sqrt(2))/2
                       0.4142135623 ); // sqrt(2)-1
  // reflections
  p = abs(p);
  p.xy -= 2.0*min(dot(vec2( k.x,k.y),p.xy),0.0)*vec2( k.x,k.y);
  p.xy -= 2.0*min(dot(vec2(-k.x,k.y),p.xy),0.0)*vec2(-k.x,k.y);
  // polygon side
  p.xy -= vec2(clamp(p.x, -k.z*r, k.z*r), r);
  vec2 d = vec2( length(p.xy)*sign(p.y), p.z-h );
  return min(max(d.x,d.y),0.0) + length(max(d,0.0));
}

float capsuleDist( vec3 p, vec3 a, vec3 b, float r )
{
        vec3 pa = p-a, ba = b-a;
        float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
        return length( pa - ba*h ) - r;
}

float roundConeDist( in vec3 p, in float r1, float r2, float h )
{
    vec2 q = vec2( length(p.xz), p.y );

    float b = (r1-r2)/h;
    float a = sqrt(1.0-b*b);
    float k = dot(q,vec2(-b,a));

    if( k < 0.0 ) return length(q) - r1;
    if( k > a*h ) return length(q-vec2(0.0,h)) - r2;

    return dot(q, vec2(a,b) ) - r1;
}

float roundConeDist2(vec3 p, vec3 a, vec3 b, float r1, float r2)
{
    // sampling independent computations (only depend on shape)
    vec3  ba = b - a;
    float l2 = dot(ba,ba);
    float rr = r1 - r2;
    float a2 = l2 - rr*rr;
    float il2 = 1.0/l2;

    // sampling dependant computations
    vec3 pa = p - a;
    float y = dot(pa,ba);
    float z = y - l2;
    float x2 = dot2( pa*l2 - ba*y );
    float y2 = y*y*l2;
    float z2 = z*z*l2;

    // single square root!
    float k = sign(rr)*rr*rr*x2;
    if( sign(z)*a2*z2 > k ) return  sqrt(x2 + z2)        *il2 - r2;
    if( sign(y)*a2*y2 < k ) return  sqrt(x2 + y2)        *il2 - r1;
                            return (sqrt(x2*a2*il2)+y*rr)*il2 - r1;
}

float triPrismDist( vec3 p, vec2 h )
{
    const float k = sqrt(3.0);
    h.x *= 0.5*k;
    p.xy /= h.x;
    p.x = abs(p.x) - 1.0;
    p.y = p.y + 1.0/k;
    if( p.x+k*p.y>0.0 ) p.xy=vec2(p.x-k*p.y,-k*p.x-p.y)/2.0;
    p.x -= clamp( p.x, -2.0, 0.0 );
    float d1 = length(p.xy)*sign(-p.y)*h.x;
    float d2 = abs(p.z)-h.y;
    return length(max(vec2(d1,d2),0.0)) + min(max(d1,d2), 0.);
}

// vertical
float cylinderVerticalDist( vec3 p, vec2 h )
{
    vec2 d = abs(vec2(length(p.xz),p.y)) - h;
    return min(max(d.x,d.y),0.0) + length(max(d,0.0));
}

// arbitrary orientation
float cylinderDist(vec3 p, vec3 a, vec3 b, float r)
{
    vec3 pa = p - a;
    vec3 ba = b - a;
    float baba = dot(ba,ba);
    float paba = dot(pa,ba);

    float x = length(pa*baba-ba*paba) - r*baba;
    float y = abs(paba-baba*0.5)-baba*0.5;
    float x2 = x*x;
    float y2 = y*y*baba;
    float d = (max(x,y)<0.0)?-min(x2,y2):(((x>0.0)?x2:0.0)+((y>0.0)?y2:0.0));
    return sign(d)*sqrt(abs(d))/baba;
}

// vertical
float coneVerticalDist( in vec3 p, in vec2 c, float h )
{
    vec2 q = h*vec2(c.x,-c.y)/c.y;
    vec2 w = vec2( length(p.xz), p.y );

        vec2 a = w - q*clamp( dot(w,q)/dot(q,q), 0.0, 1.0 );
    vec2 b = w - q*vec2( clamp( w.x/q.x, 0.0, 1.0 ), 1.0 );
    float k = sign( q.y );
    float d = min(dot( a, a ),dot(b, b));
    float s = max( k*(w.x*q.y-w.y*q.x),k*(w.y-q.y)  );
        return sqrt(d)*sign(s);
}

float cappedConeDist( in vec3 p, in float h, in float r1, in float r2 )
{
    vec2 q = vec2( length(p.xz), p.y );

    vec2 k1 = vec2(r2,h);
    vec2 k2 = vec2(r2-r1,2.0*h);
    vec2 ca = vec2(q.x-min(q.x,(q.y < 0.0)?r1:r2), abs(q.y)-h);
    vec2 cb = q - k1 + k2*clamp( dot(k1-q,k2)/dot2(k2), 0.0, 1.0 );
    float s = (cb.x < 0.0 && ca.y < 0.0) ? -1.0 : 1.0;
    return s*sqrt( min(dot2(ca),dot2(cb)) );
}

float cappedConeDist_a(vec3 p, vec3 a, vec3 b, float ra, float rb)
{
    float rba  = rb-ra;
    float baba = dot(b-a,b-a);
    float papa = dot(p-a,p-a);
    float paba = dot(p-a,b-a)/baba;

    float x = sqrt( papa - paba*paba*baba );

    float cax = max(0.0,x-((paba<0.5)?ra:rb));
    float cay = abs(paba-0.5)-0.5;

    float k = rba*rba + baba;
    float f = clamp( (rba*(x-ra)+paba*baba)/k, 0.0, 1.0 );

    float cbx = x-ra - f*rba;
    float cby = paba - f;

    float s = (cbx < 0.0 && cay < 0.0) ? -1.0 : 1.0;

    return s*sqrt( min(cax*cax + cay*cay*baba,
                       cbx*cbx + cby*cby*baba) );
}

// c is the sin/cos of the desired cone angle
float solidAngleDist(vec3 pos, vec2 c, float ra)
{
    vec2 p = vec2( length(pos.xz), pos.y );
    float l = length(p) - ra;
        float m = length(p - c*clamp(dot(p,c),0.0,ra) );
    return max(l,m*sign(c.y*p.x-c.x*p.y));
}

float octahedronDist(vec3 p, float s)
{
    p = abs(p);
    float m = p.x + p.y + p.z - s;

    // exact distance
    #if 0
    vec3 o = min(3.0*p - m, 0.0);
    o = max(6.0*p - m*2.0 - o*3.0 + (o.x+o.y+o.z), 0.0);
    return length(p - s*o/(o.x+o.y+o.z));
    #endif

    // exact distance
    #if 1
        vec3 q;
         if( 3.0*p.x < m ) q = p.xyz;
    else if( 3.0*p.y < m ) q = p.yzx;
    else if( 3.0*p.z < m ) q = p.zxy;
    else return m*0.57735027;
    float k = clamp(0.5*(q.z-q.y+s),0.0,s);
    return length(vec3(q.x,q.y-s+k,q.z-k));
    #endif

    // bound, not exact
    #if 0
        return m*0.57735027;
    #endif
}

float pyramidDist( in vec3 p, in float h )
{
    float m2 = h*h + 0.25;

    // symmetry
    p.xz = abs(p.xz);
    p.xz = (p.z>p.x) ? p.zx : p.xz;
    p.xz -= 0.5;

    // project into face plane (2D)
    vec3 q = vec3( p.z, h*p.y - 0.5*p.x, h*p.x + 0.5*p.y);

    float s = max(-q.x,0.0);
    float t = clamp( (q.y-0.5*p.z)/(m2+0.25), 0.0, 1.0 );

    float a = m2*(q.x+s)*(q.x+s) + q.y*q.y;
        float b = m2*(q.x+0.5*t)*(q.x+0.5*t) + (q.y-m2*t)*(q.y-m2*t);

    float d2 = min(q.y,-q.x*m2-q.y*0.5) > 0.0 ? 0.0 : min(a,b);

    // recover 3D and scale, and add sign
    return sqrt( (d2+q.z*q.z)/m2 ) * sign(max(q.z,-p.y));;
}

// la,lb=semi axis, h=height, ra=corner
float rhombusDist(vec3 p, float la, float lb, float h, float ra)
{
    p = abs(p);
    vec2 b = vec2(la,lb);
    float f = clamp( (ndot(b,b-2.0*p.xz))/dot(b,b), -1.0, 1.0 );
        vec2 q = vec2(length(p.xz-0.5*b*vec2(1.0-f,1.0+f))*sign(p.x*b.y+p.z*b.x-b.x*b.y)-ra, p.y-h);
    return min(max(q.x,q.y),0.0) + length(max(q,0.0));
}

float horseshoeDist( in vec3 p, in vec2 c, in float r, in float le, vec2 w )
{
    p.x = abs(p.x);
    float l = length(p.xy);
    p.xy = mat2(-c.x, c.y,
              c.y, c.x)*p.xy;
    p.xy = vec2((p.y>0.0 || p.x>0.0)?p.x:l*sign(-c.x),
                (p.x>0.0)?p.y:l );
    p.xy = vec2(p.x,abs(p.y-r))-vec2(le,0.0);

    vec2 q = vec2(length(max(p.xy,0.0)) + min(0.0,max(p.x,p.y)),p.z);
    vec2 d = abs(q) - w;
    return min(max(d.x,d.y),0.0) + length(max(d,0.0));
}

Surface unionS(Surface s1,Surface s2)
{
    if(s1.sd<=s2.sd)
        return s1;
    else
        return s2;
}

Surface intersectionS(Surface s1,Surface s2)
{
    if(s1.sd<=s2.sd)
        return s2;
    else
        return s1;
}

Surface subtractionS(Surface s1,Surface s2)
{
    Surface Ret;
    if(s1.sd<=(-s2.sd))
    {
        Ret.color=s2.color;
        Ret.Ka=s2.Ka;
        Ret.Kd=s2.Kd;
        Ret.Ks=s2.Ks;
        Ret.sd=-s2.sd;
        return Ret;
    }
    else
        return s1;
}

Surface csgObject(vec3 p)
{
    float t =iTime;
    mat4 m = rotX (30. * t) * rotY (40. * t) * rotZ (50. * t)*trans(vec3( -1.0,-0.6, -3.0)) ;
    p = opTransf(p,m);
    Surface R;
    R.sd = boxRoundDist((p), vec3(0.5),0.1);
    R.color = vec3(0.1,0.2,0.7); R.Ka=0.2; R.Kd=0.4;R.Ks=0.4;
    Surface S;
    S.sd = sphereDist(p,0.8);
    S.color = vec3(0.1,0.2,0.7); S.Ka=0.2; S.Kd=0.4;S.Ks=0.4;
    R = subtractionS(R,S);
    return R;
}
 //float ds = heartDist(p,vec3(0.0,1.0,6.0),1.0);
    // Surface Heart;
    // Heart.sd = ds; Heart.id=0;
    // Heart.color = vec3 (0.9,0.1,0.05); Heart.Ka=0.2; Heart.Kd=0.3;Heart.Ks=0.5;
    // float dp = planeDist(p,vec3 (0.0,1.0,0.0),0.0);
    // Surface Plane;
    // Plane.sd=dp;Plane.id=1;
    // Plane.color=vec3(0.75); Plane.Ka=0.2; Plane.Kd=0.4;Plane.Ks=0.4;
    // Surface d= unionS(Heart,Plane);
    // Surface CSG;
    // CSG = csgObject(p); CSG.id=2;
    // d=unionS(CSG,d);

    // Surface SphereLeft;
    // SphereLeft.sd = sphereDist(p-vec3( -3.0,1.0, 4.0),1.0);
    // SphereLeft.color = vec3(0.1,0.6,0.7); SphereLeft.Ka=0.2; SphereLeft.Kd=0.4;SphereLeft.Ks=0.4; SphereLeft.id=3;
    // d=unionS(SphereLeft,d);
    // Surface OctPrism;
    // OctPrism.sd = octogonPrismDist(p-vec3( 3.0,0.8,5.0), 0.7, 0.25); OctPrism.id=4;
    // OctPrism.color = vec3(0.4,0.2,0.7); OctPrism.Ka=0.2; OctPrism.Kd=0.7;OctPrism.Ks=0.1;
    // d=unionS(OctPrism,d);

    // Surface Sphere;
    // Sphere.id=5;
    // Sphere.sd=sphereDist(p-vec3(4.0,1.0,2.0),1.0);
    // Sphere.Ka=0.2;Sphere.Kd=0.4;Sphere.Ks=0.4;Sphere.id=5;
    // Sphere.color=vec3(0.,1.0,0.);
    // d=unionS(d,Sphere);

Surface getDist(vec3 p)
{
    Surface cyl;
    
    vec2 size = vec2(0.35, 1.2); // raio = 0.35, altura = 3.0
    vec3 position = vec3(0.0, 1.2, 4.0); // centro y = altura / 2 = 1.5
    cyl.sd = cylinderVerticalDist(p - position, size);

    //cyl.sd = cylinderVerticalDist(p - vec3(0.0, 1.0, 4.0), vec2(0.35, 1.0)); // raio = 0.5, altura = 2.0 (de -1 a 1)
    cyl.color = vec3(1.0); // cor padrão, será sobrescrita por textura
    cyl.Ka = 0.2;
    cyl.Kd = 0.4;
    cyl.Ks = 0.4;
    cyl.id = 10;  // ID novo para aplicar textura
    
    Surface ground;
    ground.sd = p.y;
    ground.color = vec3(0.2, 0.8, 0.2);
    ground.Ka = 0.1;
    ground.Kd = 0.5;
    ground.Ks = 0.1;
    ground.id = 20;

    Surface d = unionS(ground, cyl);

    // Parâmetros do sol
    float sunRadius = 0.3;
    float sunOrbitRadius = 6.0;

    // Posição animada do sol (orbita verticalmente)
    float angle = iTime * 0.2; // velocidade angular
    vec3 sunPos = vec3(sunOrbitRadius * cos(angle), 3.0 * sin(angle), 4.0);

    // Esfera do sol
    Surface sun;
    sun.sd = length(p - sunPos) - sunRadius;
    sun.color = vec3(1.0, 0.9, 0.5); // tom amarelado
    sun.Ka = 0.0;
    sun.Kd = 0.0;
    sun.Ks = 0.0;
    sun.id = 30; // use um ID exclusivo
    d = unionS(sun, d);

    return d;
}

Surface rayMarching(vec3 Cam, vec3 rd)
{
    float d0 = 0.0;
    vec3 pi = Cam + d0 * rd;
    Surface dist = getDist(pi);

    // Posição animada do sol — deve ser igual à usada em getDist
    float angle = iTime * 0.2;
    vec3 sunPos = vec3(6.0 * cos(angle), 3.0 * sin(angle), 4.0);

    float glow = 0.0;

    int i = 0;
    while ((i < MAX_STEPS) && (dist.sd > SURF_DIST) && (d0 < MAX_DIST)) {
        float distToSun = length(pi - sunPos);
        glow += exp(-distToSun * 10.0) * 0.02;

        d0 += dist.sd;
        pi = Cam + d0 * rd;
        dist = getDist(pi);
        i++;
    }

    if ((i > MAX_STEPS) || (d0 > MAX_DIST)) {
        dist.color = Sky(rd) + glow * vec3(1.0, 0.9, 0.6); // adiciona brilho ao fundo
        dist.sd = MAX_DIST;
    } else {
        dist.sd = d0;

        if (dist.id == 30) {
            dist.color = vec3(1.0, 0.9, 0.6) * 6.0; // sol visível diretamente
        }
    }

    return dist;
}

vec3 estimateNormal(vec3 p)
{
    float d= getDist(p).sd;
    float dx =  getDist(vec3(p.x+EPSILON,p.y,p.z)).sd-d;
    float dy = getDist(vec3(p.x,p.y+EPSILON,p.z)).sd-d;
    float dz = getDist(vec3(p.x-EPSILON,p.y,p.z+EPSILON)).sd-d;
    return normalize(vec3(dx,dy,dz));
}

mat3 setCamera(vec3 CamPos, vec3 LookAt)
{
   vec3 cd = normalize(LookAt-CamPos);
   vec3 cv = cross(cd,vec3(0.0,1.0,0.0));
   vec3 cu = cross(cv,cd);
   return mat3(-cv,cu,cd);
}

// https://iquilezles.org/articles/rmshadows
float calcSoftshadow(in vec3 ro, in vec3 rd, in float mint, in float tmax, int ignoreId)
{
    // bounding volume
    float tp = (tmax * 0.75 - ro.y) / rd.y;
    if (tp > 0.0) tmax = min(tmax, tp);

    float res = 1.0;
    float t = mint;

    for (int i = 0; i < 24; i++)
    {
        Surface s = getDist(ro + t * rd);

        // ignora interseções com o ID especificado (ex: o próprio sol)
        if (s.id == ignoreId) {
            t += 0.05;
            continue;
        }

        float h = s.sd;
        float sFactor = clamp(8.0 * h / t, 0.0, 1.0);
        res = min(res, sFactor);
        t += clamp(h, 0.01, 0.25);

        if (res < 0.001 || t > tmax) break;
    }

    res = clamp(res, 0.0, 1.0);
    return res * res * (3.0 - 2.0 * res); // suavização final
}


// "p" point apply texture to
// "n" normal at "p"
// "k" controls the sharpness of the blending in the
//     transitions areas.
// "s" texture sampler
vec3 boxmap( in sampler2D s, in vec3 p, in vec3 n, in float k )
{
    // project+fetch
        vec4 x = texture( s, p.yz );
        vec4 y = texture( s, p.zx );
        vec4 z = texture( s, p.xy );

    // and blend
    vec3 m = pow( abs(n), vec3(k) );
        return ((x*m.x + y*m.y + z*m.z) / (m.x + m.y + m.z)).xyz;
}

vec3 getLight(vec3 p, Surface s, vec3 Cam)
{
    if (s.sd == MAX_DIST)
        return s.color;

    // Parâmetros da luz do sol (animada)
    float angle = iTime * 0.2;
    vec3 sunPos = vec3(6.0 * cos(angle), 3.0 * sin(angle), 4.0);
    vec3 sunColor = vec3(1.0, 0.9, 0.6);

    // Segunda luz (vermelha estática)
    vec3 LightColor1 = vec3(1.0, 0.0, 0.0);
    vec3 lightPos1 = vec3(-3.0, 4.0, 4.0);

    // Vetores de direção da luz
    vec3 lightDir = normalize(sunPos - p);
    vec3 lightDir1 = normalize(lightPos1 - p);
    vec3 eye = normalize(Cam - p);

    vec3 N = estimateNormal(p);
    vec3 R = normalize(reflect(-lightDir, N));
    vec3 R1 = normalize(reflect(-lightDir1, N));

    float l = clamp(dot(N, lightDir), 0.0, 1.0);
    float l1 = clamp(dot(N, lightDir1), 0.0, 1.0);

    // Evita sombras sobre o próprio sol
    float ss = 1.0;
    float ss1 = 1.0;
    if (s.id != 30) {
        ss  = calcSoftshadow(p + 10.0 * EPSILON * lightDir,  lightDir,  0.1, 3.0, 30);
        ss1 = calcSoftshadow(p + 10.0 * EPSILON * lightDir1, lightDir1, 0.1, 3.0, -1);
    }


    // Reflexão especular
    float exp = 7.0;
    vec3 Is = vec3(0.0);
    vec3 Is1 = vec3(0.0);

    float dotRN = dot(R, eye);
    if (dotRN > 0.0)
        Is = sunColor * s.Ks * pow(dotRN, exp);

    float dotRN1 = dot(R1, eye);
    if (dotRN1 > 0.0)
        Is1 = LightColor1 * s.Ks * pow(dotRN1, exp);

    // Mapeamento de textura no cilindro
    if (s.id == 10) {
        vec3 dp = p - vec3(0.0, 1.0, 4.0);
        float theta = atan(dp.z, dp.x);
        float u = (theta + PI) / (2.0 * PI);
        float v = 1.0 - clamp(dp.y / 2.0 + 0.5, 0.0, 1.0);
        s.color = texture(iChannel0, vec2(u, v)).rgb;
    }

    // Sol visível (não iluminado, só cor intensa)
    if (s.id == 30) {
        return sunColor * 1.0;
    }

    // Composição final
    vec3 c = s.color * s.Ka +
             (s.color * sunColor * l * s.Kd) * ss +
             (s.color * LightColor1 * l1 * s.Kd) * ss1 +
             Is + Is1;

    return c;
}


    // if(s.id==3)
    // {
    //         float theta = acos(N.z);
    //         float phi = PI+ atan(N.y/N.x);
    //         float x = theta/PI;
    //         float y =phi/(2.0*PI);
    //         s.color = mix(s.color,texture(iChannel0,vec2(x,y)).xyz,0.5);
    // }

    // if(s.id==5)
    // {
    //         float theta = acos(N.z);
    //         float phi = PI+ atan(N.y/N.x);
    //         float x = theta/PI;
    //         float y =phi/(2.0*PI);
    //         float n =perlin_noise(vec2(x,y)*iResolution.xy);
    //         s.color =vec3(abs(cos(n*10.0)));

    // }

    // if(s.id==2)
    // {
    //     s.color=boxmap(iChannel0,p,N,0.5);
    // }
    


void main()
{

vec2 p = (gl_FragCoord.xy-0.5*iResolution)/iResolution.y;


vec3 rd = normalize(vec3(p,1.5));
vec3 col;
// normalized mouse coordinates


    vec2 mo = (iMouse.xy )/iResolution.xy;
    vec3 Cam= vec3 (-2.0,-1.0,1.0);
    float am = mix(-0.5*PI,0.5*PI,iMouse.z*mo.x);
    float bm = mix(-0.25*PI,0.25*PI,iMouse.z*mo.y);

    vec3 Target = vec3 (1.0,1.0,6.0);
    Cam.xz+=length(Cam-Target)*vec2(cos(am),sin(am));;
    Cam.yz+=length(Cam-Target)*vec2(cos(bm),sin(bm));;
    mat3 Ca = setCamera(Cam,Target);
    rd=Ca*rd;

Surface d = rayMarching(Cam,rd);
vec3 po= Cam+d.sd*rd;
vec3 l = getLight(po,d,Cam);

// col = l;
if (d.sd < MAX_DIST) {
    col = l; // cor da superfície com iluminação
} else {
    col =  vec3(0.4, 0.7, 1.0); // céu
    // if (rd.y < -0.01) {
    //     col = vec3(0.2, 0.6, 0.2); // grama
    // } else {
    //     col = vec3(0.4, 0.7, 1.0); // céu
    // }
}


// gamma
col = pow( col, vec3(0.4545) );
C = vec4( col, 1.0 );

}

