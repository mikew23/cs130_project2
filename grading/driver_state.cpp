#include "driver_state.h"
#include <cstring>

driver_state::driver_state()
{
}

driver_state::~driver_state()
{
    delete [] image_color;
    delete [] image_depth;
}

// This function should allocate and initialize the arrays that store color and
// depth.  This is not done during the constructor since the width and height
// are not known when this class is constructed.
void initialize_render(driver_state& state, int width, int height)
{
    state.image_width=width;
    state.image_height=height;
    //state.image_color=0;
    state.image_depth= new float[width*height];
    std::cout<<"TODO: allocate and initialize state.image_color and state.image_depth."<<std::endl;
    state.image_color = new pixel[width* height];
    for(int i = 0; i < width * height; i++){
    	state.image_color[i] = make_pixel(0,0,0);
	state.image_depth[i] = 1;
    }
}

// This function will be called to render the data that has been stored in this class.
// Valid values of type are:
//   render_type::triangle - Each group of three vertices corresponds to a triangle.
//   render_type::indexed -  Each group of three indices in index_data corresponds
//                           to a triangle.  These numbers are indices into vertex_data.
//   render_type::fan -      The vertices are to be interpreted as a triangle fan.
//   render_type::strip -    The vertices are to be interpreted as a triangle strip.
void render(driver_state& state, render_type type)
{
    std::cout<<"TODO: implement rendering."<<std::endl;
    switch(type){
    	case render_type::triangle:{
		data_geometry g[3];
		data_vertex v[3];
		for (int i = 0; i< state.num_vertices*state.floats_per_vertex; i+=3*state.floats_per_vertex){
        		v[0].data = &state.vertex_data[i];
        		v[1].data = &state.vertex_data[i+state.floats_per_vertex*1];
        		v[2].data = &state.vertex_data[i+state.floats_per_vertex*2];
      		}

		for (int i = 0; i < 3; i++){
        		state.vertex_shader(v[i], g[i], state.uniform_data);
                }
       		rasterize_triangle(state, g[0], g[1],g[2]);
    		
		//delete [] g;
		//delete [] v;
		break;
	}
   	case render_type::indexed:{
		break;
	}
	case render_type::fan:{
		break;
	}
	case render_type::strip:{
		break;
	}
	default:{
		break;
	} 
         			
    }
}


// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry& v0,
    const data_geometry& v1, const data_geometry& v2,int face)
{
    if(face==6)
    {
        rasterize_triangle(state, v0, v1, v2);
        return;
    }
    std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
    clip_triangle(state, v0, v1, v2,face+1);
}

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.
void rasterize_triangle(driver_state& state, const data_geometry& v0,
    const data_geometry& v1, const data_geometry& v2)
{
    std::cout<<"TODO: implement rasterization"<<std::endl;
/*
  vec4 positiona = v0.gl_Position/v0.gl_Position[3];
  vec4 positionb = v1.gl_Position/v1.gl_Position[3];
  vec4 positionc = v2.gl_Position/v2.gl_Position[3];
  float w0 = v0.gl_Position[3];
  float w1 = v1.gl_Position[3];
  float w2 = v2.gl_Position[3];

  float x0 = v0.gl_Position[0] / w0;
  float x1 = v1.gl_Position[0] / w1;
  float x2 = v2.gl_Position[0] / w2;

  float y0 = v0.gl_Position[1] / w0;
  float y1 = v1.gl_Position[1] / w1;
  float y2 = v2.gl_Position[1] / w2;

  float z0 = v0.gl_Position[2] / w0;
  float z1 = v1.gl_Position[2] / w1;
  float z2 = v2.gl_Position[2] / w2;
  vec3 Z = {z0, z1, z2};

  float Ai = (state.image_width/2) * x0 + (state.image_width/2 - 0.5);
  float Aj = (state.image_height/2)* y0 + (state.image_height/2 - 0.5);
  float Bi = (state.image_width/2) * x1 + (state.image_width/2 - 0.5);
  float Bj = (state.image_height/2)* y1 + (state.image_height/2 - 0.5);
  float Ci = (state.image_width/2) * x2 + (state.image_width/2 - 0.5);
  float Cj = (state.image_height/2)* y2 + (state.image_height/2 - 0.5);


*/   
    float Ai = (state.image_width/2) * (v0.gl_Position[0]/v0.gl_Position[3]+1);
    float Aj = (state.image_height/2)* (v0.gl_Position[1]/v0.gl_Position[3]+1);
    float Bi = (state.image_width/2) * (v1.gl_Position[0]/v1.gl_Position[3]+1);
    float Bj = (state.image_height/2)* (v1.gl_Position[1]/v1.gl_Position[3]+1);
    float Ci = (state.image_width/2) * (v2.gl_Position[0]/v2.gl_Position[3]+1);
    float Cj = (state.image_height/2)* (v2.gl_Position[1]/v2.gl_Position[3]+1);
/*
    float Ai = ((state.image_width* (v0.gl_Position[0]/v0.gl_Position[3]))+ state.image_width -1)/2;
    float Aj = ((state.image_height/2)* (v0.gl_Position[1]/v0.gl_Position[3])+(state.image_height/2) -0.5);
    float Bi = ((state.image_width/2) * (v1.gl_Position[0]/v1.gl_Position[3])+(state.image_width/2) -0.5);
    float Bj = ((state.image_height/2)* (v1.gl_Position[1]/v1.gl_Position[3])+(state.image_height/2) -0.5);
    float Ci = ((state.image_width/2) * (v2.gl_Position[0]/v2.gl_Position[3])+(state.image_width/2) -0.5);
    float Cj = ((state.image_height/2)* (v2.gl_Position[1]/v2.gl_Position[3])+(state.image_height/2) -0.5);
   */ 
    int x_min = std::min(std::min(Ai,Bi),Ci);
    int x_max = std::max(std::max(Ai,Bi),Ci);
    int y_min = std::min(std::min(Aj,Bj),Cj);
    int y_max = std::max(std::max(Aj,Bj),Cj);

    //keep x, y min and max within the limit
    x_min = (x_min < 0? 0 : x_min);
    x_max = (x_max > state.image_width? state.image_width -1 : x_max);
    y_min = (y_min < 0? 0 : y_min);
    y_max = (y_max > state.image_height? state.image_height -1 : y_max);

   // float alpha_per, beta_per, gamma_per;
    int image_index;
      vec2 A = {Ai,Aj};
  vec2 B = {Bi,Bj};
  vec2 C = {Ci,Cj};
    float area_ABC =  (((Bi * Cj) - (Ci * Bj)) - ((Ai * Cj) - (Ci * Aj)) + ((Ai * Bj) - (Bi * Aj))) * 0.5;
    //float totalArea = (C[0] - A[0]) * (B[1] - A[1]) - (C[1] - A[1]) * (B[0] - A[0]);
    //loop through all the pixels
    for (int i = x_min; i < x_max; i++){
	for (int j = y_min; j < y_max ; j++){
     		float a = 0.5 * (((Bi * Cj) - (Ci * Bj)) - ((i * Cj) - (Ci * j)) + ((i * Bj) - (Bi * j)));
      		float b = 0.5 * (((i * Cj) - (Ci * j)) - ((Ai * Cj) - (Ci * Aj)) + ((Ai * j) - (i * Aj)));
    	  	float g = 0.5 * (((Bi * j) - (i * Bj)) - ((Ai * j) - (i * Aj)) + ((Ai * Bj) - (Bi * Aj)));
   	   	a = a / area_ABC;
		b= b / area_ABC;
		g= g / area_ABC;
	/*
		vec2 P = {(float)i,(float)j};
		float a = ((P[0] - B[0]) * (C[1] - B[1]) - (P[1] - B[1]) * (C[0] - B[0]))/totalArea;
		float b = ((P[0] - C[0]) * (A[1] - C[1]) - (P[1] - C[1]) * (A[0] - C[0]))/totalArea;
 		float g = ((P[0] - A[0]) * (B[1] - A[1]) - (P[1] - A[1]) * (B[0] - A[0]))/totalArea;
		*/
		image_index = i + j * state.image_width;
		if (a >= 0 && b >= 0 && g<= 1){
        		int index = state.image_width * j + i;
        		state.image_color[index] = make_pixel(255,255,255);
	/*		data_fragment frag;
			frag.data = new float[MAX_FLOATS_PER_VERTEX];
                	data_output output;
			float depth1 = alpha * v0.gl_Position[2]/v0.gl_Position[3] + beta * v1.gl_Position[2]/v1.gl_Position[3] + gamma * v2.gl_Position[2]/v2.gl_Position[3];
      			if (state.image_depth[image_index] > depth1) {
                    for (int k = 0; k < state.floats_per_vertex; k++) {
                        float k_gour;
                        switch (state.interp_rules[k]) {
                            case interp_type::flat:
                                frag.data[k] = v0.data[k];
                            break;

                            case interp_type::smooth:
                                k_gour = (alpha/v0.gl_Position[3] + beta/v1.gl_Position[3] + gamma/v2.gl_Position[3]);

                                alpha_per = alpha/k_gour/(v0.gl_Position[3]);
                                beta_per = beta/k_gour/(v1.gl_Position[3]);
                                gamma_per = gamma/k_gour/(v2.gl_Position[3]);

                                frag.data[k] = alpha_per * v0.data[k] + beta_per * v1.data[k] + gamma_per * v2.data[k];
                            break;

                            case interp_type::noperspective:
                                frag.data[k] = alpha * v0.data[k] + beta * v1.data[k] + gamma * v2.data[k];

                            break;

                            default:
                            break;
                        }
                    }
                    state.fragment_shader(frag, output, state.uniform_data);
                    output.output_color = output.output_color * 255;
                    state.image_color[image_index] = make_pixel(output.output_color[0], output.output_color[1], output.output_color[2]);
                    state.image_depth[image_index] = depth1;
                }*/
		}
		float alpha, beta, gamma;
/*
      if (a >= 0 && b >= 0 && g>= 0)
      {
        data_output out_data;
        data_fragment fragment;
        fragment.data = new float[state.floats_per_vertex];

        float depth = (a * Z[0]) + (b * Z[1]) + (g* Z[2]);

        if(state.image_depth[image_index] > depth)
        {
        for(int i = 0; i < state.floats_per_vertex; i++)
        {

          if(state.interp_rules[i] == interp_type::flat)
          {
            fragment.data[i] = v0.data[i];
          }
          if(state.interp_rules[i] == interp_type::smooth)
          {
		float k = (a/ w0) + (b/ w1) + (g/ w2);
		alpha = a / w0 / k;
            beta = b/ w1 / k;
            gamma = g/ w2 / k;
            fragment.data[i] = (alpha*v0.data[i])
            + (beta*v1.data[i])
            + (gamma*v2.data[i]);
          }
          if(state.interp_rules[i] == interp_type::noperspective)
          {
            fragment.data[i] = a*v0.data[i]
                             + b*v1.data[i]
                             + g*v2.data[i];
          }
          else if(state.interp_rules[i] == interp_type::invalid)
          {
            std::cout << "Error" << std::endl;
          }
        }
      }

        state.fragment_shader(fragment, out_data, state.uniform_data);
        out_data.output_color = out_data.output_color*255;
        int red   = out_data.output_color[0];
        int green = out_data.output_color[1];
        int blue  = out_data.output_color[2];

        float z = a*positiona[2]+b*positionb[2]+g*positionc[2];
        int index = (state.image_width * j) + i;
	if(state.image_depth[index] > z)
        {state.image_depth[index] = z;
          state.image_color[index] = make_pixel(red,green,blue);
	}
    	}*/
	}
    }
}


