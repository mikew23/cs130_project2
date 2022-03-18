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
		for (int i = 0; i < state.num_vertices; i += 3) { 
                data_geometry* tri_array = new data_geometry[3];
                for (int j = 0; j < 3; j++) { 
                    data_vertex ver;
                    ver.data = new float[MAX_FLOATS_PER_VERTEX];
                    tri_array[j].data = new float[MAX_FLOATS_PER_VERTEX];
                    for (int k = 0; k < state.floats_per_vertex; k++) {
                        ver.data[k] = state.vertex_data[k + state.floats_per_vertex*(i+j)];
                        tri_array[j].data[k] = ver.data[k];
                    }
                    state.vertex_shader((const data_vertex)ver, tri_array[j], state.uniform_data);
                }
		/*data_geometry g[3];
		data_vertex v[3];
		for (int i = 0; i< state.num_vertices*state.floats_per_vertex; i+=3*state.floats_per_vertex){
        		v[0].data = &state.vertex_data[i];
        		v[1].data = &state.vertex_data[i+state.floats_per_vertex*1];
        		v[2].data = &state.vertex_data[i+state.floats_per_vertex*2];
      		}

		for (int i = 0; i < 3; i++){
        		state.vertex_shader(v[i], g[i], state.uniform_data);
                }*/
       		rasterize_triangle(state, tri_array[0], tri_array[1],tri_array[2]);
    	}	
		//delete [] g;
		//delete [] v;
		break;
	}
   	case render_type::indexed:{
		/*data_geometr*y dg_array = new data_geometry[3];
		data_vertex dv_array[3];
		for (int i = 0; i < state.num_triangles * 3; i += 3) {
            		for (int j = 0; j < 3; j++) {
                		dv_array[j].data = state.vertex_data[state.index_data[i + j] * state.floats_per_vertex];
		                dg_array[j].data = dv_array[j].data;
               			state.vertex_shader(dv_array[j], dg_array[j], state.uniform_data);
 				}
            	clip_triangle(state, dg_array[0], dg_array[1], dg_array[2], 0);
        	}*/
		break;
	}
	case render_type::fan:{
		/* int flag;
		data_geometry* dg_tri = new data_geometry[3];
    		data_vertex dv[3];
        	for (int i = 0; i < state.num_vertices; i++) {
            		for (int j = 0; j < 3; j++) {
                		flag = i + j;
                		if (j == 0) flag = 0;
                		dv[j].data = &state.vertex_data[flag * state.floats_per_vertex];
                		dg_tri[j].data = dv[j].data;
                		state.vertex_shader(dv[j], dg_tri[j], state.uniform_data);
             	   		
           		 }
            		clip_triangle(state, dg_temp[0], dg_temp[1], dg_temp[2], 0);
       		 }*/
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
    //std::cout<<"TODO: implement rasterization"<<std::endl;
/*    float Ai = (state.image_width/2.0) * (in[0]->gl_Position[0]/in[0]->gl_Position[3]) + (state.image_width/2.0) - (0.5);
    float Aj = (state.image_height/2.0) * (in[0]->gl_Position[1]/in[0]->gl_Position[3]) + (state.image_height/2.0) - (0.5);

    float Bi = (state.image_width/2.0) * (in[1]->gl_Position[0]/in[1]->gl_Position[3]) + (state.image_width/2.0) - (0.5);
    float Bj = (state.image_height/2.0) * (in[1]->gl_Position[1]/in[1]->gl_Position[3]) + (state.image_height/2.0) - (0.5);

    float Ci = (state.image_width/2.0) * (in[2]->gl_Position[0]/in[2]->gl_Position[3]) + (state.image_width/2.0) - (0.5);
    float Cj = (state.image_height/2.0) * (in[2]->gl_Position[1]/in[2]->gl_Position[3]) + (state.image_height/2.0) - (0.5); 

    float Ai = (state.image_width/2) * (v0.gl_Position[0]/v0.gl_Position[3]+1);
    float Aj = (state.image_height/2)* (v0.gl_Position[1]/v0.gl_Position[3]+1);
    float Bi = (state.image_width/2) * (v1.gl_Position[0]/v1.gl_Position[3]+1);
    float Bj = (state.image_height/2)* (v1.gl_Position[1]/v1.gl_Position[3]+1);
    float Ci = (state.image_width/2) * (v2.gl_Position[0]/v2.gl_Position[3]+1);
    float Cj = (state.image_height/2)* (v2.gl_Position[1]/v2.gl_Position[3]+1);
*/
    float Ai = ((state.image_width* (v0.gl_Position[0]/v0.gl_Position[3]))+ state.image_width -1)/2;
    float Aj = ((state.image_height/2)* (v0.gl_Position[1]/v0.gl_Position[3])+(state.image_height/2) -0.5);
    float Bi = ((state.image_width/2) * (v1.gl_Position[0]/v1.gl_Position[3])+(state.image_width/2) -0.5);
    float Bj = ((state.image_height/2)* (v1.gl_Position[1]/v1.gl_Position[3])+(state.image_height/2) -0.5);
    float Ci = ((state.image_width/2) * (v2.gl_Position[0]/v2.gl_Position[3])+(state.image_width/2) -0.5);
    float Cj = ((state.image_height/2)* (v2.gl_Position[1]/v2.gl_Position[3])+(state.image_height/2) -0.5);
float a_temp, b_temp, g_temp;   
    int x_min = std::min(std::min(Ai,Bi),Ci);
    int x_max = std::max(std::max(Ai,Bi),Ci);
    int y_min = std::min(std::min(Aj,Bj),Cj);
    int y_max = std::max(std::max(Aj,Bj),Cj);
    float area_ABC =  (((Bi * Cj) - (Ci * Bj)) - ((Ai * Cj) - (Ci * Aj)) + ((Ai * Bj) - (Bi * Aj))) * 0.5;
    //keep x, y min and max within the limit
    x_min = (x_min < 0? 0 : x_min);
    y_min = (y_min < 0? 0 : y_min);
    x_max = (x_max > state.image_width  ? state.image_width  - 1 : x_max);
    y_max = (y_max > state.image_height ? state.image_height - 1 : y_max);
     for (int i = x_min; i < x_max; i++){
	for (int j = y_min; j < y_max ; j++){
     		float a = 0.5 * (((Bi * Cj) - (Ci * Bj)) - ((i * Cj) - (Ci * j)) + ((i * Bj) - (Bi * j)));
      		float b = 0.5 * (((i * Cj) - (Ci * j)) - ((Ai * Cj) - (Ci * Aj)) + ((Ai * j) - (i * Aj)));
    	  	float g = 0.5 * (((Bi * j) - (i * Bj)) - ((Ai * j) - (i * Aj)) + ((Ai * Bj) - (Bi * Aj)));
   	   	a = a / area_ABC;
		b = b / area_ABC;
		g = g / area_ABC;
/*		if (a >= 0 && b >= 0 && g<= 1){
        		int index = state.image_width * j + i;
        		state.image_color[index] = make_pixel(255,255,255);
      		}*/
		if (a>= 0 && b >= 0 && g >= 0) {
		 int index = state.image_width * j + i;
		data_fragment frag; //{data};
                frag.data = new float[MAX_FLOATS_PER_VERTEX];
                data_output output;
		float depth1 = a * v0.gl_Position[2]/ v0.gl_Position[3] + b* v1.gl_Position[2]/v1.gl_Position[3] + g * v2.gl_Position[2]/v2.gl_Position[3];

                if (state.image_depth[index] > depth1) {
                    for (int k = 0; k < state.floats_per_vertex; k++) {
                        float k_temp;
                        switch (state.interp_rules[k]) {
                            case interp_type::flat:
                                frag.data[k] = v0.data[k];
                            break;

                            case interp_type::smooth:
                                k_temp = (a/v0.gl_Position[3] + b/v1.gl_Position[3] + g/v2.gl_Position[3]);
                                a_temp = a/k_temp/(v0.gl_Position[3]);
                                b_temp = b/k_temp/(v1.gl_Position[3]);
                                g_temp = g/k_temp/(v2.gl_Position[3]);

                                frag.data[k] = a_temp * v0.data[k] + b_temp * v1.data[k] + g_temp * v2.data[k];
                            break;

                            case interp_type::noperspective:
                                frag.data[k] = a * v0.data[k] + b * v1.data[k] + g * v2.data[k];

                            break;

                            default:
                            break;
			}
                    state.fragment_shader(frag, output, state.uniform_data);
                    output.output_color = output.output_color * 255;
                    state.image_color[index] = make_pixel(output.output_color[0], output.output_color[1], output.output_color[2]);
                    state.image_depth[index] = depth1;
			}
		}
	}
	}
    }
}


