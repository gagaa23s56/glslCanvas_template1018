// Author:CMH
// Title:BreathingGlow
#ifdef GL_ES
precision mediump float;
#endif

#define PI 3.14159
#define M_SQRT_2 1.41421356

uniform vec2 u_resolution;
uniform vec2 u_mouse;
uniform float u_time;


vec3 hash( vec3 p ) // replace this by something better
    {
        p = vec3( dot(p,vec3(127.1,311.7, 74.7)),
                dot(p,vec3(269.5,183.3,246.1)),
                dot(p,vec3(113.5,271.9,124.6)));

        return -1.0 + 2.0*fract(sin(p)*43758.5453123);
    }

//Gradient Noise
vec2 hash2( vec2 x )        //亂數範圍 [-1,1]
{
    const vec2 k = vec2( 0.3183099, 0.3678794 );
    x = x*k + k.yx;
    return -1.0 + 2.0*fract( 16.0 * k*fract( x.x*x.y*(x.x+x.y)) );
}
float gnoise( in vec2 p )   //亂數範圍 [-1,1]
{
    vec2 i = floor( p );
    vec2 f = fract( p );

    vec2 u = f*f*(3.0-2.0*f);

    return mix( mix( dot( hash2( i + vec2(0.0,0.0) ), f - vec2(0.0,0.0) ), 
        dot( hash2( i + vec2(1.0,0.0) ), f - vec2(1.0,0.0) ), u.x),
        mix( dot( hash2( i + vec2(0.0,1.0) ), f - vec2(0.0,1.0) ), 
        dot( hash2( i + vec2(1.0,1.0) ), f - vec2(1.0,1.0) ), u.x), u.y);
}
float noise( in vec3 p )
    {
        vec3 i = floor( p );
        vec3 f = fract( p );
        
        vec3 u = f*f*(3.0-2.0*f);

        return mix( mix( mix( dot( hash( i + vec3(0.0,0.0,0.0) ), f - vec3(0.0,0.0,0.0) ), 
                        dot( hash( i + vec3(1.0,0.0,0.0) ), f - vec3(1.0,0.0,0.0) ), u.x),
                    mix( dot( hash( i + vec3(0.0,1.0,0.0) ), f - vec3(0.0,1.0,0.0) ), 
                        dot( hash( i + vec3(1.0,1.0,0.0) ), f - vec3(1.0,1.0,0.0) ), u.x), u.y),
                mix( mix( dot( hash( i + vec3(0.0,0.0,1.0) ), f - vec3(0.0,0.0,1.0) ), 
                        dot( hash( i + vec3(1.0,0.0,1.0) ), f - vec3(1.0,0.0,1.0) ), u.x),
                    mix( dot( hash( i + vec3(0.0,1.0,1.0) ), f - vec3(0.0,1.0,1.0) ), 
                        dot( hash( i + vec3(1.0,1.0,1.0) ), f - vec3(1.0,1.0,1.0) ), u.x), u.y), u.z );
    }
float fbm(in vec2 uv)       //亂數範圍 [-1,1]
{
    float f;                //fbm - fractal noise (4 octaves)
    mat2 m = mat2( 1.6,  1.2, -1.2,  1.6 );
    f   = 0.5000*gnoise( uv ); uv = m*uv;		  
    f += 0.2500*gnoise( uv ); uv = m*uv;
    f += 0.1250*gnoise( uv ); uv = m*uv;
    f += 0.0625*gnoise( uv ); uv = m*uv;
    return f;
}

float random (vec2 st) {
        return fract(sin(dot(st.xy, vec2(12.9898,78.233)))*43758.5453123);
}

float glow(float d, float str, float thickness){
    return thickness / pow(d, str);
}

float square(vec2 P, float size){
        return abs(max(abs(P.x), abs(P.y)) - size/(2.0*M_SQRT_2));
    }

mat2 rotate2d(float _angle){
        return mat2(cos(_angle),-sin(_angle), sin(_angle),cos(_angle));
}

float disc(vec2 P, float size){
        return length(P) - size/2.;
}

float heart(vec2 P, float size){
        float x = M_SQRT_2/2. * (P.x - P.y);
        float y = M_SQRT_2/2. * (P.x + P.y);
        float r1 = max(abs(x),abs(y))-size/3.5;
        float r2 = length(M_SQRT_2/2.*vec2(1.,-1.)-size/3.5);
        
        float r3 = length(M_SQRT_2/2.*vec2(-1.,-1.)-size/3.5);
        
        return min(min(r1,r2),r3);
    }


float sdFish(float i, vec2 p, float a) {
        float ds, c = cos(a), s = sin(a);
        p *= 1.3*mat2(c,s,-s,c); // Rotate and rescale
        p.x *= .97 + (.04+.2*p.y)*cos(i+9.*u_time);  
        // Swiming ondulation (+rotate in Z axes)
        
        ds = min(length(p-vec2(.8,0))-.45, length(p-vec2(-.14,0))-.12);   
        // Distance to fish

        p.y = abs(p.y)+.13;
        return max(min(length(p),length(p-vec2(.56,0)))-.3,-ds)*.05;
}

float mouseEffect(vec2 uv, vec2 mouse, float size)
{
    float dist=length(uv-mouse);
    return 1.-smoothstep(size*10., size, dist);  //size
    //return pow(dist, 0.5);
}

float circle(vec2 uv, float radius){
    float dist = length(uv);
    float circle_dist = abs(dist-radius);                                //光環大小
    return circle_dist;
}

void main() {
    vec2 uv = gl_FragCoord.xy/u_resolution.xy;
    uv= uv*2.0-1.0;
    uv.x *= u_resolution.x/u_resolution.y;
    
    //grid repetition
    // vec2 uvs=uv*2.;
    // vec2 ipos = floor(uvs);  // get the integer coords
    // vec2 fpos = fract(uvs);  // get the fractional coords
    // uv= fpos*2.0-1.0;
    // uv*= rotate2d(random(ipos)*PI+u_time*0.5);
    vec3 color = vec3(0.245,0.087,0.137);
    
    vec2 mouse=u_mouse/u_resolution.xy;
    mouse=mouse*2.0-1.0;
    mouse.x*= u_resolution.x/u_resolution.y;
    
    float interact=1.-mouseEffect(uv,mouse,0.08);
    
    //陰晴圓缺
        float theta=2.0*PI*u_time/8.0;
        vec2 point=vec2(sin(theta), cos(theta));
        float dir= dot(point, (uv))+0.55;
        float theta2=2.0*PI*u_time/12.0;
        vec2 point2=vec2(sin(theta2), cos(theta2));
        float dir2 = (dot(point2, (uv))+0.55)*.2+.8;
    
    int total = 12;
    
    for(int index=0;index<12;++index)
    {
//         float noise_position= smoothstep(0.2, 0.7, uv.x+0.420);
//         // float noise_position= interact;
//         float noise_scale=0.328*noise_position;
//         float noise_freq=2.084;
//         //float circle_dist=circle(uv, 0.480+noise_scale
//             //*noise( vec3(uv*noise_freq, float(index)+u_time*0.4)) );

//         uv *= rotate2d(-0.004);
//         // float circle_dist=square(uv, 0.480+noise_scale
//         //     *noise( vec3(uv*noise_freq, float(index)+u_time*0.4)) );
        
        float noise_position= interact;
        float radius_noise=noise(vec3(4.892*uv,float(index)+u_time*0.388))*0.280*noise_position;
        float radius=0.572+radius_noise;
        float circle_dist = circle(uv, radius);   


        //定義圓環
        float dist = length(uv)+0.176*fbm(uv*1.520*dir + vec2(index));
        // float circle_dist = abs(dist-0.512);								//光環大小

        //動態呼吸
        float breathing=sin(u_time*2.0* PI /7.424)*0.5+0.5;						//option1
        //float breathing=(exp(sin(u_time/2.0*pi)) - 0.36787944)*0.42545906412; 			//option2 正確
        //float strength =(0.2*breathing*dir+0.180);			//[0.2~0.3]			//光暈強度加上動態時間營造呼吸感
        float strength =(0.2*breathing+-0.076);			//[0.2~0.3]			//光暈強度加上動態時間營造呼吸感
        float thickness=(0.1*breathing+0.084);			//[0.1~0.2]			//光環厚度 營造呼吸感
        float glow_circle = glow(circle_dist, strength, thickness)/float(total)*1.8;



        color += glow_circle; //*vec3(1.000,0.348,0.006)
    }
    
    //亂數作用雲霧
    float fog= (fbm(0.4*uv+vec2(-0.920*u_time, -0.02*u_time))*0.6+0.1)*.15;
    
    
    gl_FragColor = vec4(vec3(color)*vec3(0.975,0.511,0.455)*dir2-fog,1.0);
    // gl_FragColor = vec4(color,1.0);
}

