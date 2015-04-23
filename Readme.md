Nemo is an individual-based forward-time simulation program created and maintained by Fred Guillaume. The latest release is available for download here: http://nemo2.sourceforge.net/

This repository contains source files for Nemo which have been modified with the help of Fred to perform two main new functions: allow breeding to occur within a specified window, i.e. not only within a given patch, and to allow a range-wide dispersal kernel along with a specified connectivity matrix to be provided in the input rather than a full patch number by patch number matrix of dispersal probabilities. Details follow below.

The current version compiles with Nemo version 2.3.4. I have not tested whether they compile compatibly with v2.2.0 (current release; new release of v2.3.0 is expected soon). Preliminary tests have been done to ensure the code works as desired.

Nemo also requires that the GNU Scientific Library be installed in order to compile and run. Brian O'Meara's lab page has a simple how-to for installing the GSL: http://www.brianomeara.info/tutorials/brownie/gsl

## Documentation

#### Dispersal

Two new parameters in the .init file allow for a constant dispersal kernel across patches to be input.

1: "dispersal_kernel" accepts a single array of sorted (or unsorted, though this will result in longer run times) dispersal probabilities (forward migration) from the natal patch. The patches that these probabilities correspond to will be defined in the connectivity matrix. This specified dispersal kernel is recycled across all patches. If one wishes to have this vary across patches, Fred's original version can be re-implemented, or there are comments in LCEdisperse.cc where these modifications should be made.

2: "dispersal_connectivity_matrix" accepts a matrix with the number of rows corresponding to the number of patches and the number of columns corresponding to the length of the dispersal kernel. These values will be the patch IDs (integer numbers) that correspond to the probabilities given in the dispersal kernel. Row 1 contains the patches that patch 1 will send migrants to, row 2 for patch 2, etc.

#### Breeding Window

Four new parameters in the .init file allow for use of a breeding window and a couple related features. Similar to the dispersal kernel, the breeding window will let fathers be chosen from the specified patches with the specified probabilities. This can be set to happen every time an offspring is created, or set to only occur if the patch containing the female has no males present. In the latter case, if there is no male in the focal patch, and there is a female, and still no male is found within the breeding window, the female can be told to self.

1: "self_if_alone" is a boolean parameter, which if specified in the .init file will tell females to self if they are unable to find a male to mate with.

2: "always_breed_window" is a boolean parameter, which if specified in the .init file will ensure that a breeding window is used every single time an offspring is created. Otherwise, the breeding window will only be used if a male cannot be found in the focal patch containing the female.

3: "breeding_kernel" takes the same format as "dispersal_kernel" above, except that these probabilities are the forward probabilities that any other nearby patch would have of contributing a father to the focal patch. Within the code, these probabilities are then normalized based on the number of males present at any given time in the potential patches.

4: "breeding_connectivity_matrix" takes the same format as "dispersal_connectivity_matrix" above, and contains the patch IDs of all patches within the breeding window, corresponding to the probabilities provided in "breeding_kernel".

#### Deleterious Mutations

One small modification has been made to the gamma distribution for deleterious mutations, where 30% of the time, a mutation is lethal. Do not include this version of ttdeletmutations_bitstring.cc if you would not like this condition to be imposed.
