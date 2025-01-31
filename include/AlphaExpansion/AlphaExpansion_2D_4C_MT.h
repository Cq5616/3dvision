// Alpha-Expansion solver for discrete multi-label optimization on grids
// Written by Lenka Saidlova at the Czech Technical University in Prague

// This software is in the public domain. Where that dedication is not
// recognized, you are granted a perpetual, irrevocable license to copy
// and modify this file however you want.

#ifndef AlphaExpansion_2D_4C_MT_H_
#define AlphaExpansion_2D_4C_MT_H_

#include <algorithm>
#include <GridCut/GridGraph_2D_4C_MT.h>

//0:  [ 1, 0]
//1:  [ 0, 1]
#define ALPHAEXPANSION_NEIGHBORS 2
#define ALPHAEXPANSION_INFINITY 1000000

template<
	typename type_label, // Type used to represent labels
	typename type_cost,  // Type used to represent data and smoothness costs
	typename type_energy // Type used to represent resulting energy
>
class AlphaExpansion_2D_4C_MT
{

public:
	// Function for representing the smoothness term.
	typedef type_energy (*SmoothCostFn)(int, int, int, int);
	SmoothCostFn smooth_fn;

	// Constructors of the algorithm.
	// See "README.txt" how to fill data and smoothness costs
	AlphaExpansion_2D_4C_MT(int width, int height, int nLabels, type_cost *data, type_cost **smooth, int num_threads, int block_size);
	AlphaExpansion_2D_4C_MT(int width, int height, int nLabels, type_cost *data, type_cost *smooth, int num_threads, int block_size);
	AlphaExpansion_2D_4C_MT(int width, int height, int nLabels, type_cost *data, SmoothCostFn smooth_fn, int num_threads, int block_size);
	~AlphaExpansion_2D_4C_MT(void);
	
	// Sets the all pixels to the given label
	void set_labels(type_label label);

	// Sets the new labeling.
	// The array has to have width*height elements, one value for each pixel.
	void set_labeling(type_label* labeling);

	// Runs the minimization. 
	// The algorithm iterate over the labels in a fixed order, from 0 to k - 1
	// and it stops when no further improvement is possible.
	void perform();

	// Runs the minimization. 
	// The algorithm iterate over the labels in a fixed order, from 0 to k - 1
	// and it stops after the number of max_cycles.
	void perform(int max_cycles);

	// Runs the minimization. 
	// The algorithm iterate over the labels in a random order
	// and it stops when no further improvement is possible.
	void perform_random();

	// Runs the minimization. 
	// The algorithm iterate over the labels in a random order
	// and it stops after the number of max_cycles.
	void perform_random(int max_cycles);

	// Returns the resulting energy.
	type_energy get_energy();

	// Returns the final labeling.
	// The array has width*height elements, one value for each pixel.
	type_label* get_labeling();

	// Returns the final label for the given coordinates
	type_label get_label(int x, int y);

	// Returns the final label for the given pixel index
	type_label get_label(int pix);

private:
	int width;
	int height;
	int nLabels; 
	long int nPixels; 
	type_cost **smooth_cost;
	type_cost *data_cost;
	type_label *labeling;
	bool smooth_vary;
	bool smooth_array;
	typedef GridGraph_2D_4C_MT<type_cost,type_cost,type_cost> Grid;
	Grid* grid;
	int num_threads;
	int block_size;
	
	type_cost get_data_cost(int pix,int lab) const;
	type_cost get_smooth_cost(int pix,int lab1,int lab2) const;
	void perform(int max_cycles, bool random);
	void init_common(int width, int height, int nLabels, type_cost *data, int num_threads, int block_size);
	type_energy perform_cycle(bool random);
	void create_grid_array(type_label alpha_label);
	void create_grid_fn(type_label alpha_label);
	void perform_label(type_label alpha_label);
	type_energy get_data_energy();
	type_energy get_smooth_energy();
	type_energy get_smooth_energy_vary();
	type_energy get_smooth_energy_fn();	
	void add_tlink(int pix, type_cost to_source, type_cost to_sink, type_cost *source, type_cost *sink);
	void add_h_nlink(int pix, int nPix, type_cost A, type_cost B, type_cost C, type_cost D, 
		type_cost *source, type_cost *sink, type_cost *left, type_cost *top, type_cost *right, type_cost *bottom);
	void add_v_nlink(int pix, int nPix, type_cost A, type_cost B, type_cost C, type_cost D, 
		type_cost *source, type_cost *sink, type_cost *left, type_cost *top, type_cost *right, type_cost *bottom);
};

template<typename type_label, typename type_cost, typename type_energy>
type_cost AlphaExpansion_2D_4C_MT<type_label, type_cost, type_energy>::get_data_cost(int pix,int lab) const
{
	return data_cost[(pix)*nLabels + lab];
}

template<typename type_label, typename type_cost, typename type_energy>
type_cost AlphaExpansion_2D_4C_MT<type_label, type_cost, type_energy>::get_smooth_cost(int pix,int lab1,int lab2) const
{
	return smooth_cost[pix][(lab2)+(lab1)*nLabels]; 
}


template<typename type_label, typename type_cost, typename type_energy>
void AlphaExpansion_2D_4C_MT<type_label, type_cost, type_energy>::init_common(int width, int height, int nLabels, type_cost *data, int num_threads, int block_size)
{
	this->width = width;
	this->height = height;
	this->nLabels = nLabels;
	this->nPixels = width*height;
	this->data_cost = data;
	this->num_threads = num_threads;
	this->block_size = block_size;
	
	labeling = new type_label[nPixels];
	std::fill(labeling, labeling + nPixels, 0);

}

template<typename type_label, typename type_cost, typename type_energy>
AlphaExpansion_2D_4C_MT<type_label, type_cost, type_energy>::AlphaExpansion_2D_4C_MT(int width, int height, int nLabels, type_cost *data, type_cost **smooth, int num_threads, int block_size)
{

	init_common(width, height, nLabels, data, num_threads, block_size);
	
	smooth_vary = true;
	smooth_array = true;
	smooth_cost = smooth;
	
}

template<typename type_label, typename type_cost, typename type_energy>
AlphaExpansion_2D_4C_MT<type_label, type_cost, type_energy>::AlphaExpansion_2D_4C_MT(int width, int height, int nLabels, type_cost *data, type_cost *smooth, int num_threads, int block_size)
{
	init_common(width, height, nLabels, data, num_threads, block_size);
		
	smooth_vary = false;
	smooth_array = true;
	smooth_cost = new type_cost*[1];
	smooth_cost[0] = smooth;
	
}

template<typename type_label, typename type_cost, typename type_energy>
AlphaExpansion_2D_4C_MT<type_label, type_cost, type_energy>::AlphaExpansion_2D_4C_MT(int width, int height, int nLabels, type_cost *data, SmoothCostFn smoothFn, int num_threads, int block_size)
{
	init_common(width, height, nLabels, data, num_threads, block_size);
		
	smooth_vary = false;
	smooth_array = false;
	smooth_cost = NULL;
	this->smooth_fn = smoothFn;
	
}

template<typename type_label, typename type_cost, typename type_energy>
void AlphaExpansion_2D_4C_MT<type_label, type_cost, type_energy>::set_labeling(type_label* new_labeling)
{

	this->labeling = new_labeling;
}

template<typename type_label, typename type_cost, typename type_energy>
void AlphaExpansion_2D_4C_MT<type_label, type_cost, type_energy>::set_labels(type_label label)
{

	std::fill(labeling, labeling + nPixels, label);
}

template<typename type_label, typename type_cost, typename type_energy>
type_energy AlphaExpansion_2D_4C_MT<type_label, type_cost, type_energy>::get_data_energy()
{

	type_energy energy = (type_energy) 0;

	for ( int i = 0; i < nPixels; i++ ){
		energy += get_data_cost(i, labeling[i]);
	}

	return energy;
}

template<typename type_label, typename type_cost, typename type_energy>
type_energy AlphaExpansion_2D_4C_MT<type_label, type_cost, type_energy>::get_smooth_energy_vary()
{

	type_energy energy = (type_energy) 0;
	int pix;

	for ( int y = 0; y < height; y++ ){
		for ( int x = 0; x < width - 1; x++ ){
			pix = x+y*width;
			energy += get_smooth_cost(ALPHAEXPANSION_NEIGHBORS*pix, labeling[pix], labeling[pix+1]);
		}
	}

	for ( int y = 0; y < height - 1; y++ ){
		for ( int x = 0; x < width; x++ ){
			pix = x+y*width;
			energy += get_smooth_cost(ALPHAEXPANSION_NEIGHBORS*pix+1, labeling[pix], labeling[pix+width]);
		}
	}

	return energy;
}

template<typename type_label, typename type_cost, typename type_energy>
type_energy AlphaExpansion_2D_4C_MT<type_label, type_cost, type_energy>::get_smooth_energy_fn()
{

	type_energy energy = (type_energy) 0;
	int pix, nPix;

	for ( int y = 0; y < height; y++ ){
		for ( int x = 0; x < width - 1; x++ ){
			pix = x+y*width;
			nPix = pix+1;
			energy += smooth_fn(pix, nPix, labeling[pix], labeling[nPix]);
		}
	}

	for ( int y = 0; y < height - 1; y++ ){
		for ( int x = 0; x < width; x++ ){
			pix = x+y*width;
			nPix = pix+width;
			energy += smooth_fn(pix, nPix, labeling[pix], labeling[nPix]);
		}
	}

	return energy;
}

template<typename type_label, typename type_cost, typename type_energy>
type_energy AlphaExpansion_2D_4C_MT<type_label, type_cost, type_energy>::get_smooth_energy()
{

	type_energy energy = (type_energy) 0;
	int pix;

	for ( int y = 0; y < height; y++ ){
		for ( int x = 0; x < width - 1; x++ ){
			pix = x+y*width;
			energy += get_smooth_cost(0, labeling[pix], labeling[pix+1]);
		}
	}

	for ( int y = 0; y < height - 1; y++ ){
		for ( int x = 0; x < width; x++ ){
			pix = x+y*width;
			energy += get_smooth_cost(0, labeling[pix], labeling[pix+width]);
		}
	}

	return energy;
}

template<typename type_label, typename type_cost, typename type_energy>
type_energy AlphaExpansion_2D_4C_MT<type_label, type_cost, type_energy>::get_energy()
{

	if(smooth_array){
		if(smooth_vary){
			return get_data_energy() + get_smooth_energy_vary();
		}
		return get_data_energy() + get_smooth_energy();
	}
	return get_data_energy() + get_smooth_energy_fn();
}

template<typename type_label, typename type_cost, typename type_energy>
void AlphaExpansion_2D_4C_MT<type_label, type_cost, type_energy>::perform()
{

	perform(ALPHAEXPANSION_INFINITY, false);
}

template<typename type_label, typename type_cost, typename type_energy>
void AlphaExpansion_2D_4C_MT<type_label, type_cost, type_energy>::perform(int max_cycles)
{

	perform(max_cycles, false);
}

template<typename type_label, typename type_cost, typename type_energy>
void AlphaExpansion_2D_4C_MT<type_label, type_cost, type_energy>::perform_random(int max_cycles)
{

	perform(max_cycles, true);
}

template<typename type_label, typename type_cost, typename type_energy>
void AlphaExpansion_2D_4C_MT<type_label, type_cost, type_energy>::perform_random(){

	perform(ALPHAEXPANSION_INFINITY, true);
}

template<typename type_label, typename type_cost, typename type_energy>
void AlphaExpansion_2D_4C_MT<type_label, type_cost, type_energy>::perform(int max_cycles, bool random)
{

	int cycle = 1;
    type_energy new_energy = get_energy();
	type_energy old_energy = -1;

	while ((old_energy < 0 || old_energy > new_energy)  && cycle <= max_cycles){
        old_energy = new_energy;
        new_energy = perform_cycle(random);
		cycle++;   
    }

}

template<typename type_label, typename type_cost, typename type_energy>
type_energy AlphaExpansion_2D_4C_MT<type_label, type_cost, type_energy>::perform_cycle(bool random)
{

	int* order = new int[nLabels];
	for (int i = 0; i < nLabels; i++){ 
		order[i] = i;
	}
	
	if(random){	
		std::random_shuffle(order, order + nLabels);	
	}

	for (int i = 0;  i < nLabels;  i++ ){
		perform_label(order[i]);
	}
       
    return get_energy();
}

template<typename type_label, typename type_cost, typename type_energy>
void AlphaExpansion_2D_4C_MT<type_label, type_cost, type_energy>::perform_label(type_label alpha_label)
{

	grid = new Grid(width,height,num_threads,block_size);
	if(smooth_array){
		create_grid_array(alpha_label);
	}else{
		create_grid_fn(alpha_label);
	}
	grid->compute_maxflow();
	
	int pix = 0;
	for (int y = 0; y < height; y++){
		for (int x = 0; x < width; x++){
            
			if ( labeling[pix] != alpha_label ){
                if (!grid->get_segment(grid->node_id(x,y))){
                    labeling[pix] = alpha_label;
				}
            }
			pix++;
        }
	}	
	
	delete grid;
	
}

template<typename type_label, typename type_cost, typename type_energy>
void AlphaExpansion_2D_4C_MT<type_label, type_cost, type_energy>::create_grid_array(type_label alpha_label)
{

	int pix, nPix;
	int* indices = new int[ALPHAEXPANSION_NEIGHBORS];
	std::fill(indices, indices + ALPHAEXPANSION_NEIGHBORS, 0);
	
	type_cost *source = new type_cost[nPixels];
	type_cost *sink = new type_cost[nPixels];
	type_cost *left = new type_cost[nPixels];
	type_cost *top = new type_cost[nPixels];
	type_cost *right = new type_cost[nPixels];
	type_cost *bottom = new type_cost[nPixels];
	std::fill(source, source + nPixels, 0);
	std::fill(sink, sink + nPixels, 0);
	std::fill(left, left + nPixels, 0);
	std::fill(top, top + nPixels, 0);
	std::fill(right, right + nPixels, 0);
	std::fill(bottom, bottom + nPixels, 0);
	
	for (int y = 0; y < height; y++){
		for (int x = 0; x < width; x++){

			pix = y*width + x;
			if(smooth_vary){
				for (int i = 0; i < ALPHAEXPANSION_NEIGHBORS; i++){
					indices[i] = ALPHAEXPANSION_NEIGHBORS*pix + i;
				}
			}
			
			if(labeling[pix] != alpha_label){

				add_tlink(pix, get_data_cost(pix, labeling[pix]), get_data_cost(pix, alpha_label), source, sink);	

				if (x < width - 1){				
					nPix = pix + 1;
					if(labeling[nPix] != alpha_label){
						add_h_nlink(pix, nPix, get_smooth_cost(indices[0], alpha_label, alpha_label), get_smooth_cost(indices[0], alpha_label, labeling[nPix]), 
							get_smooth_cost(indices[0], labeling[pix], alpha_label), get_smooth_cost(indices[0], labeling[pix], labeling[nPix]), 
							source, sink, left, top, right, bottom);
					}else{
						add_tlink(pix, get_smooth_cost(indices[0], labeling[pix], alpha_label), get_smooth_cost(indices[0], alpha_label, labeling[nPix]), source, sink);
					}
				}
				if(y < height - 1){
					nPix = pix + width;
					if(labeling[nPix] != alpha_label ){
						add_v_nlink(pix, nPix, get_smooth_cost(indices[1], alpha_label, alpha_label), get_smooth_cost(indices[1], alpha_label, labeling[nPix]), 
							get_smooth_cost(indices[1], labeling[pix], alpha_label), get_smooth_cost(indices[1], labeling[pix], labeling[nPix]), 
							source, sink, left, top, right, bottom);
					}else{
						add_tlink(pix, get_smooth_cost(indices[1], labeling[pix], alpha_label), get_smooth_cost(indices[1], alpha_label, labeling[nPix]), source, sink);
					}
				}
			}else{
				if (x < width - 1){				
					nPix = pix + 1;
					if(labeling[nPix] != alpha_label){
						add_tlink(nPix, get_smooth_cost(indices[0], alpha_label, labeling[nPix]), get_smooth_cost(indices[0], alpha_label, alpha_label), source, sink);
					}
				}
				if(y < height - 1){
					nPix = pix + width;
					if(labeling[nPix] != alpha_label){
						add_tlink(nPix, get_smooth_cost(indices[1], alpha_label, labeling[nPix]), get_smooth_cost(indices[1], alpha_label, alpha_label), source, sink);
					}
				}
			}			
		}
	}	

	grid->set_caps(source, sink, left, right, top, bottom);

	delete [] source;
	delete [] sink;
	delete [] left;
	delete [] top;
	delete [] right;
	delete [] bottom;
	delete [] indices;

}

template<typename type_label, typename type_cost, typename type_energy>
void AlphaExpansion_2D_4C_MT<type_label, type_cost, type_energy>::create_grid_fn(type_label alpha_label)
{

	int pix, nPix;
	type_cost *source = new type_cost[nPixels];
	type_cost *sink = new type_cost[nPixels];
	type_cost *left = new type_cost[nPixels];
	type_cost *top = new type_cost[nPixels];
	type_cost *right = new type_cost[nPixels];
	type_cost *bottom = new type_cost[nPixels];
	std::fill(source, source + nPixels, 0);
	std::fill(sink, sink + nPixels, 0);
	std::fill(left, left + nPixels, 0);
	std::fill(top, top + nPixels, 0);
	std::fill(right, right + nPixels, 0);
	std::fill(bottom, bottom + nPixels, 0);

	for (int y = 0; y < height; y++){
		for (int x = 0; x < width; x++){
	
			pix = y*width + x;
						
			if(labeling[pix] != alpha_label){

				add_tlink(pix, get_data_cost(pix, labeling[pix]), get_data_cost(pix, alpha_label), source, sink);	

				if (x < width - 1){				
					nPix = pix + 1;
					if(labeling[nPix] != alpha_label){
						add_h_nlink(pix, nPix, smooth_fn(pix, nPix, alpha_label, alpha_label), smooth_fn(pix, nPix, alpha_label, labeling[nPix]), 
							smooth_fn(pix, nPix, labeling[pix], alpha_label), smooth_fn(pix, nPix, labeling[pix], labeling[nPix]), 
							source, sink, left, top, right, bottom);
					}else{
						add_tlink(pix, smooth_fn(pix, nPix, labeling[pix], alpha_label), smooth_fn(pix, nPix, alpha_label, labeling[nPix]), source, sink);
					}
				}
				if(y < height - 1){
					nPix = pix + width;
					if(labeling[nPix] != alpha_label ){
						add_v_nlink(pix, nPix, smooth_fn(pix, nPix, alpha_label, alpha_label), smooth_fn(pix, nPix, alpha_label, labeling[nPix]), 
							smooth_fn(pix, nPix, labeling[pix], alpha_label), smooth_fn(pix, nPix, labeling[pix], labeling[nPix]), 
							source, sink, left, top, right, bottom);
					}else{
						add_tlink(pix, smooth_fn(pix, nPix, labeling[pix], alpha_label), smooth_fn(pix, nPix, alpha_label, labeling[nPix]), source, sink);
					}
				}
			}else{
				if (x < width - 1){				
					nPix = pix + 1;
					if(labeling[nPix] != alpha_label){
						add_tlink(nPix, smooth_fn(pix, nPix, alpha_label, labeling[nPix]), smooth_fn(pix, nPix, alpha_label, alpha_label), source, sink);
					}
				}
				if(y < height - 1){
					nPix = pix + width;
					if(labeling[nPix] != alpha_label){
						add_tlink(nPix, smooth_fn(pix, nPix, alpha_label, labeling[nPix]), smooth_fn(pix, nPix, alpha_label, alpha_label), source, sink);
					}
				}
			}			
		}
	}	
	
	grid->set_caps(source, sink, left, right, top, bottom);

	delete [] source;
	delete [] sink;
	delete [] left;
	delete [] top;
	delete [] right;
	delete [] bottom;
}

template<typename type_label, typename type_cost, typename type_energy>
void AlphaExpansion_2D_4C_MT<type_label, type_cost, type_energy>::add_h_nlink(int pix, int nPix, type_cost A, type_cost B, type_cost C, type_cost D, 
																			type_cost *source, type_cost *sink, type_cost *left, type_cost *top, type_cost *right, type_cost *bottom)
{
	
	if ( A+D > C+B) {
		type_cost delta = A+D-C-B;
        type_cost subtrA = delta/3;

        A = A-subtrA;
        C = C+subtrA;
        B = B+(delta-subtrA*2);
	}
	add_tlink(pix, D, A, source, sink);
    B -= A; 
	C -= D;
	
    if (B < 0){
		
		add_tlink(pix, -B, 0, source, sink);
		add_tlink(nPix, 0, -B, source, sink);

		left[nPix] += B+C;

    }else if (C < 0){
        
		add_tlink(pix, 0, -C, source, sink);
		add_tlink(nPix, -C, 0, source, sink);

		right[pix] += B+C;
		
    }else{
		
		right[pix] += B;
		left[nPix] += C;
    }
	
}

template<typename type_label, typename type_cost, typename type_energy>
void AlphaExpansion_2D_4C_MT<type_label, type_cost, type_energy>::add_v_nlink(int pix, int nPix, type_cost A, type_cost B, type_cost C, type_cost D, 
																			type_cost *source, type_cost *sink, type_cost *left, type_cost *top, type_cost *right, type_cost *bottom)
{

	if ( A+D > C+B) {
		type_cost delta = A+D-C-B;
        type_cost subtrA = delta/3;

        A = A-subtrA;
        C = C+subtrA;
        B = B+(delta-subtrA*2);
	}
	add_tlink(pix, D, A, source, sink);
    B -= A; 
	C -= D;
	
    if (B < 0){
		
		add_tlink(pix, -B, 0, source, sink);
		add_tlink(nPix, 0, -B, source, sink);

		top[nPix] += B+C;

    }else if (C < 0){
        
		add_tlink(pix, 0, -C, source, sink);
		add_tlink(nPix, -C, 0, source, sink);

		bottom[pix] += B+C;
	
    }else{
		
		bottom[pix] += B;
		top[nPix] += C;
    }
	
}

template<typename type_label, typename type_cost, typename type_energy>
void AlphaExpansion_2D_4C_MT<type_label, type_cost, type_energy>::add_tlink(int pix, type_cost to_source, type_cost to_sink, type_cost *source, type_cost *sink)
{

	source[pix] += to_source;
	sink[pix] += to_sink;
}

template<typename type_label, typename type_cost, typename type_energy>
type_label * AlphaExpansion_2D_4C_MT<type_label, type_cost, type_energy>::get_labeling(void)
{

	return labeling;
}

template<typename type_label, typename type_cost, typename type_energy>
type_label AlphaExpansion_2D_4C_MT<type_label, type_cost, type_energy>::get_label(int x, int y)
{

	return labeling[y*width + x];
}

template<typename type_label, typename type_cost, typename type_energy>
type_label AlphaExpansion_2D_4C_MT<type_label, type_cost, type_energy>::get_label(int pix)
{

	return labeling[pix];
}

template<typename type_label, typename type_cost, typename type_energy>
AlphaExpansion_2D_4C_MT<type_label, type_cost, type_energy>::~AlphaExpansion_2D_4C_MT(void)
{
	delete [] labeling;
	delete [] data_cost;
	if(smooth_array)
		delete [] smooth_cost;

}

#undef ALPHAEXPANSION_NEIGHBORS
#undef ALPHAEXPANSION_INFINITY

#endif
