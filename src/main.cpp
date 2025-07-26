//
// Created by Alexander on 6/14/2025.
//

// OpenGL-related includes
#include <glad/glad.h>
#include <glfw/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/type_ptr.hpp>

// stl includes
#include <iostream>
#include <vector>
#include <fstream>
#include <chrono>
// throw a SIGINT to stop iterations

// local includes
#include "shader.h"
#include "math_util.h"

void print_timestamp(const std::chrono::time_point<std::chrono::steady_clock>& start) {
    auto end = std::chrono::steady_clock::now();
    auto duration = end - start;
    auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds> (duration);

    printf("Elapsed time::PRINT_TIMESTAMP: %lld ms\n", elapsed_ms.count());
}

// TODO: TMP VARIABLES, ADD A CONSTEXPR EVALUATION, MIGRATE THESE TO MATH_UTIL.H
constexpr int CURL_ROWS = 98;
constexpr int CURL_COLS = 398;

void framebuffer_size_callback(GLFWwindow *window, int width, int height)
{
    glViewport(0,0,width,height);
}

void processInput(GLFWwindow *window)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
    {
        glfwSetWindowShouldClose(window, true);
    }
}

bool fileExists(const char* fileName)
{
    std::ifstream test(fileName);
    if (test) {
        return true;
    } else {
        return false;
    }
}

// VBO vertices: fullscreen quad
constexpr float vertices[] = {
        // Positions                    // texCoords
        -1.0f,  1.0f, 0.0f, 0.0f, 1.0f,
        1.0f,  1.0f, 0.0f, 1.0f, 1.0f,
        1.0f, -1.0f, 0.0f, 1.0f, 0.0f,
        -1.0f, -1.0f, 0.0f, 0.0f, 0.0f
};

// EBO to collapse 6 vertices into 4
constexpr unsigned int indices[] = {
        0, 1, 2,  // First triangle
        2, 3, 0   // Second triangle
};


int main() {
    // glfw and shader init
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);

    GLFWwindow  *window = glfwCreateWindow(WIDTH,HEIGHT, "grid", nullptr,nullptr);

    if (window == nullptr)
    {
        std::cout << "Failed to create new window" << std::endl;
        glfwTerminate();
        exit(1);
    }

    glfwMakeContextCurrent(window);

    if (!gladLoadGLLoader((GLADloadproc) glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

    glViewport(0,0,WIDTH,HEIGHT);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

    // for tracking fps
    double lastTime = glfwGetTime();
    int frameCount = 0;

    try {
        assert(fileExists(VERT_SHADER_PATH));
        assert(fileExists(FRAG_SHADER_PATH));
    } catch (const std::exception& e) {
        std::perror(e.what());
        exit(1);
    }

    Shader shaderProgram = Shader(VERT_SHADER_PATH, FRAG_SHADER_PATH);

    // end glfw and shader init
    ThreeD F(Ny, MatrixXd(Nx, NL));

    for (MatrixXd& mat : F) {
        constexpr static double coeff{0.01};
        mat.fill(1);
        mat.setConstant(true);
        MatrixXd rand(Nx, NL);
        rand.setRandom();
        rand = rand * coeff;

        mat = mat + rand;

        // set rightward velocity
        for (int i = 0; i < mat.rows(); i++) {
            mat(i, 3) = RIGHT_VEL;
        }
    }

    MatrixXb cylinder = MatrixXb(Ny, Nx);
    cylinder.setConstant(true);
    cylinder.fill(false);
    int bndry_sz{0};

    for (int y = 0; y < Ny; y++) {
        for (int x = 0; x < Nx; x++) {
            // only dynamic_cast incurs overhead
            // static_cast should still be 3 extra asm instructions
            if (squared_distance(static_cast<float>(x), static_cast<float>(y), Nx / 4., Ny / 2.) < CYLINDER_RADIUS_SQUARED) {
                cylinder(y,x) = true;
                bndry_sz++;
            }
        }
    }

    // continue openGL setup
    unsigned int VAO;
    glGenVertexArrays(1, &VAO);
    glBindVertexArray(VAO);

    unsigned int VBO;
    glGenBuffers(1, &VBO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);

    unsigned int EBO;
    glGenBuffers(1, &EBO);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);

    // VAO ATTRIBUTES
    // Position attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    // texCoord (location = 1)
    glVertexAttribPointer(1, 2,GL_FLOAT, GL_FALSE, 5 * sizeof(float), (void*)(3*sizeof(float)));
    glEnableVertexAttribArray(1);

    // use shader
    GLuint texture_ID;
    glGenTextures(1, &texture_ID);                // gen texture id
    glActiveTexture(GL_TEXTURE0);                     // Activate texture unit 0
    glBindTexture(GL_TEXTURE_2D, texture_ID);   // bind texture

    // shader parameters
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    // Allocate texture memory
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, CURL_COLS, CURL_ROWS, 0, GL_RED, GL_FLOAT, nullptr);

    //  setup shader + uniform
    glUseProgram(shaderProgram.ID);
    // fetch uniform
    // dataTexture
    GLint loc = glGetUniformLocation(shaderProgram.ID, "dataTexture");
    glUniform1i(loc, 0); // set sampler2D to use GL_TEXTURE0
    // minValue
    GLint loc2 = glGetUniformLocation(shaderProgram.ID, "minValue");
    glUniform1f(loc2, 0.0f);
    // maxValue
    GLint loc3 = glGetUniformLocation(shaderProgram.ID, "maxValue");
    glUniform1f(loc3, 0.0f);

    // main loop
    // TODO: set this to infinite loop with signal mask for SIGINT from ctrl+c
    while (!glfwWindowShouldClose(window)) {
        // Zou-He absorb boundary conditions
        // smooth out fluid pressure at boundary by introducing curl?
        auto start = std::chrono::steady_clock::now();
        // TODO: rewrite F dimensions order
        /*
         * The last position represents timeseries window/snapshot of values at a specific timestep, this should be what's controlled by the vector container
         * We need to make F as (NL, Ny, Nx), NOT (Ny, Nx, NL)
         * Same case for Feq
         * This allows more intuitive code and more "native" Eigen function calls which will greatly reduce overhead
         * ::--Changes to Make--::
         * F init, Feq init
         * Boundary conditions: can access matrixXd by idx
         * Streaming step: Self-explanatory
         * Apply boundary (F[cylinder, :])
         * Feq collision
         */

        for (MatrixXd &mat: F) {
            mat(Ny - 1, 6) = mat(Ny - 2, 6);
            mat(Ny - 1, 7) = mat(Ny - 2, 7);
            mat(Ny - 1, 8) = mat(Ny - 2, 8);

            mat(0, 2) = mat(1, 2);
            mat(0, 3) = mat(1, 3);
            mat(0, 4) = mat(1, 4);
        }

        for (int i = 0; i < NL; i++) {
            int cx = cxs[i];
            int cy = cys[i];

            // streaming step: roll
            roll(F, cx, 1, i);
            roll(F, cy, 0, i);
            // this is probably most prone to optimization
        }

        print_timestamp(start);

        // apply bndry mask
        MatrixXd bndryF = apply_boundary(F, cylinder, bndry_sz);
        // reorder lattice points by column
        matrix_reorder(bndryF, {0, 5, 6, 7, 8, 1, 2, 3, 4});
        MatrixXd rho = sum_axis_two(F);

        // pass by value or reference in sum_axis_two? currently pass by ref but two LOC
        ThreeD tmp_x = element_prod(F, cxs);
        MatrixXd tmp_sum_x = sum_axis_two(tmp_x);
        MatrixXd ux = element_div(tmp_sum_x, rho);

        ThreeD tmp_y = element_prod(F, cys);
        MatrixXd tmp_sum_y = sum_axis_two(tmp_y);
        MatrixXd uy = element_div(tmp_sum_y, rho);

        // apply boundary
        apply_bndry_to_F(F, cylinder, bndryF);
        apply_boundary_to_vel(ux, cylinder);
        apply_boundary_to_vel(uy, cylinder);

        // collision
        ThreeD Feq = ThreeD(Ny, MatrixXd(Nx, NL));
        for (MatrixXd &mat: Feq) {
            mat.setZero();
        }

        for (int i = 0; i < NL; i++) {
            int cx = cxs[i], cy = cys[i];
            double w = weights[i];

            // collision update
//            IFUCKEDUP(ux);
//            IFUCKEDUP(uy);
            MatrixXd weighted_square_matrix = (cx * ux + cy * uy).array().square();
            MatrixXd ux_sq = ux.array().square();
            MatrixXd uy_sq = uy.array().square();
            MatrixXd sum_sq = ux_sq + uy_sq;
            MatrixXd collision_val = (rho * w).array() *
                                     (1. + 3. * (cx * ux + cy * uy).array() + 9. * (weighted_square_matrix).array() / 2. +
                                      (sum_sq).array() / 2.);

//            IFUCKEDUP(collision_val);
            for (int i_feq = 0; i_feq < Feq.size(); i_feq++) {
                MatrixXd& mat = Feq[i_feq];
                for (int j = 0; j < mat.rows(); j++) {
                    mat.row(j).setConstant(collision_val(i_feq,j));
                }
            }
        }
        // F = F  +-(1/tau) * (F-Feq)
        for (int i = 0; i < F.size(); i++) {
            auto &Fi = F[i];
            auto &Feqi = Feq[i];
            Fi = Fi.array() + (-1. / tau) * (Fi - Feqi).array();
        }

        // calculate curl and plot
        MatrixXd curl = get_curl(ux, uy);
        // TODO: add shader to window by passing it into a VAO with a full-screen quad
        // TODO: Add OpenGL
        static const int cols = static_cast<int>(curl.cols());
        static const int rows = static_cast<int>(curl.rows());
        // TMP ASSERTION UNTIL WE FIX CONSTEXPR
        assert(cols == CURL_COLS && rows == CURL_ROWS);
        Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> float_matrix = curl.cast<float>();

        auto end = std::chrono::steady_clock::now();
        auto duration = end - start;
        auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds> (duration);

        printf("Elapsed time: %lld ms\n", elapsed_ms.count());

        // fps
        double currentTime = glfwGetTime();
        frameCount++;
        if (currentTime - lastTime >= 1.0) {
            double fps = static_cast<double>(frameCount) / (currentTime - lastTime);
            printf("Current FPS: %.2f\n", fps);

            frameCount = 0;
            lastTime = currentTime;
        }

        // opengl updates
        glUseProgram(shaderProgram.ID);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, texture_ID);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, cols, rows, 0, GL_RED, GL_FLOAT, float_matrix.data());

        // draw
        glClearColor(0.2f,0.3f,0.3f,1.0f);
        glBindVertexArray(VAO);
        glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, nullptr);
        // clear + swap
        glfwSwapBuffers(window);
        glClear(GL_COLOR_BUFFER_BIT);
        glfwPollEvents();
    }
    return 0;
}