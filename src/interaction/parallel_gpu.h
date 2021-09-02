#include <gsl/gsl_vector.h>
#include "src/utils/print_utils.h"
#include "src/wake/wake.h"
#include "src/interaction/interaction.h"
#include "src/interaction/interaction_utils.h"
#include "src/utils/timing_utils.h"

namespace pawan{
    class __parallel_gpu : public __interaction{

    private:

        //! Interact
        /*
         * Compute interaction between particles of a single wake object
         * \param	W	Wake object pointer
         */
        virtual void interact(__wake *W);                           //there is no need to call derived class methods as virtual (remove them everywhere?)

        //! Interact
        /*
         * Compute interaction between particles of two wake objects
         * \param	W1	Wake 1 object pointer
         * \param	W2	Wake 2 object pointer
         */
        virtual void interact(__wake *W1, __wake *W2);

    public:

        //! Constructor
        /*
         * Creates single wake interaction object
         * \param	W	Wake object pointer
         */
        __parallel_gpu(__wake *W);

        //! Constructor
        /*
         * Creates empty interaction object with two wakes
         * \param	W1	Wake object pointer
         * \param	W2	Wake object pointer
         */
        __parallel_gpu(__wake *W1, __wake *W2);

        //! Destructor
        /*
         * Deletes particles
         */
        ~__parallel_gpu() = default;
    };
}
#endif
