#ifndef MAIN_UTILS_HPP
#define MAIN_UTILS_HPP

#include "mesh/Mesh_config.h"
#include <iostream>

// =============================================================================
// run_mesh_inspector
//
// Interactive menu that lets the user inspect any of the three staggered
// meshes (P, U, V) after generation. Call this at the end of main() if
// post-generation inspection is desired.
// =============================================================================
inline void run_mesh_inspector(const MeshConfig& config,
                                const GridData&   pMesh,
                                const GridData&   uMesh,
                                const GridData&   vMesh)
{
    char option;
    bool running = true;

    while (running) {
        std::cout << "\n=========================================" << std::endl;
        std::cout << "      MESH INSPECTION SELECTOR           " << std::endl;
        std::cout << "=========================================" << std::endl;
        std::cout << " [P] Inspect P Mesh (Pressure/Centers)"  << std::endl;
        std::cout << " [U] Inspect U Mesh (Velocity X)"        << std::endl;
        std::cout << " [V] Inspect V Mesh (Velocity Y)"        << std::endl;
        std::cout << " [Q] Quit"                               << std::endl;
        std::cout << "-----------------------------------------" << std::endl;
        std::cout << " Select an option: ";
        std::cin >> option;
        option = toupper(option);

        switch (option) {
            case 'P': config.inspect_mesh(pMesh); break;
            case 'U': config.inspect_mesh(uMesh); break;
            case 'V': config.inspect_mesh(vMesh); break;
            case 'Q':
                std::cout << "Exiting inspector..." << std::endl;
                running = false;
                break;
            default:
                std::cout << "Invalid option. Please try again." << std::endl;
        }
    }
}

#endif // MAIN_UTILS_HPP
