
#include "sphere_collision.hpp"

#include <random>

#ifdef SCENE_SPHERE_COLLISION

using namespace vcl;







void scene_model::frame_draw(std::map<std::string,GLuint>& , scene_structure& scene, gui_structure& )
{
    float dt = 0.02f * timer.scale;
    timer.update();

    set_gui();

    create_new_particle();
    compute_time_step(dt);

    display_particles(scene);
    draw(borders, scene.camera);

}

void scene_model::compute_time_step(float dt)
{
    // Set forces
    const size_t N = particles.size();
    for(size_t k=0; k<N; ++k)
        particles[k].f = vec3(0,-9.81f,0);


    // Integrate position and speed of particles through time
    for(size_t k=0; k<N; ++k) {
        particle_structure& particle = particles[k];
        vec3& v = particle.v;
        vec3& p = particle.p;
        vec3 const& f = particle.f;

        v = (1-0.9f*dt) * v + dt * f; // gravity + friction force
        p = p + dt * v;
    }

    // Cube normals
    std::vector<vec3> normals = {
        vec3( 0,  0, -1),  // Blue, z
        vec3( 0, -1,  0),  // Green, y
        vec3(-1,  0,  0),  // Red, x
        vec3( 0,  0,  1),  // -Blue, -z
        vec3( 0,  1,  0),  // -Green, -y
        vec3( 1,  0,  0)   // -Red, -x
    };

    for (size_t i = 0; i < N; ++i) {
        particle_structure& particle = particles[i];  // Current particle

        // Collisions between spheres
        for (size_t j = 0; j < N; ++j) {
            if (j == i) { continue; }
            particle_structure& particle2 = particles[j];

            float d = norm(particle.p - particle2.p);
            if (d < 2 * particle.r) {  // Collision
                vec3 dir = (particle2.p - particle.p) / d;
                // if (dot(dir, ))
                particle2.p += dir * (2*particle.r - d);  // Update position
                particle2.v += dir * norm(particle.v) * 0.8;  // Update velocity
            }
        }

        // Collisions with cube
        for (size_t k = 0; k < 6; ++k) {
            vec3 n = normals[k];
            vec3 a = particle.p - vec3(
                    particle.p[0] * abs(n[0]),
                    particle.p[1] * abs(n[1]),
                    particle.p[2] * abs(n[2])
            ) - n; // Projection
            float detection = dot(particle.p - a, n);

            if (detection <= particle.r) {
                float d = norm(particle.p - a);
                particle.p += (particle.r - d) * n;  // Reset position
                particle.v += n * norm(particle.v) * 0.8;  // Reset velocity
            }
        }
    }
}


void scene_model::create_new_particle()
{
    // Emission of new particle if needed
    timer.periodic_event_time_step = gui_scene.time_interval_new_sphere;
    const bool is_new_particle = timer.event;
    static const std::vector<vec3> color_lut = {
        {1,0,0},{0,1,0},{0,0,1},{1,1,0},{1,0,1},{0,1,1}
    };

    if( is_new_particle && gui_scene.add_sphere)
    {
        particle_structure new_particle;

        new_particle.r = 0.08f;
        new_particle.c = color_lut[int(rand_interval()*color_lut.size())];

        // Initial position
        new_particle.p = vec3(0,0,0);

        // Initial speed
        const float theta = rand_interval(0, 2*3.14f);
        new_particle.v = vec3( 2*std::cos(theta), 5.0f, 2*std::sin(theta));

        particles.push_back(new_particle);

    }
}
void scene_model::display_particles(scene_structure& scene)
{
    const size_t N = particles.size();
    for(size_t k=0; k<N; ++k)
    {
        const particle_structure& part = particles[k];

        sphere.uniform.transform.translation = part.p;
        sphere.uniform.transform.scaling = part.r;
        sphere.uniform.color = part.c;
        draw(sphere, scene.camera);
    }
}




void scene_model::setup_data(std::map<std::string,GLuint>& shaders, scene_structure& , gui_structure& )
{
    sphere = mesh_drawable( mesh_primitive_sphere(1.0f));
    sphere.shader = shaders["mesh"];

    std::vector<vec3> borders_segments = {{-1,-1,-1},{1,-1,-1}, {1,-1,-1},{1,1,-1}, {1,1,-1},{-1,1,-1}, {-1,1,-1},{-1,-1,-1},
                                          {-1,-1,1} ,{1,-1,1},  {1,-1,1}, {1,1,1},  {1,1,1}, {-1,1,1},  {-1,1,1}, {-1,-1,1},
                                          {-1,-1,-1},{-1,-1,1}, {1,-1,-1},{1,-1,1}, {1,1,-1},{1,1,1},   {-1,1,-1},{-1,1,1}};
    borders = segments_gpu(borders_segments);
    borders.uniform.color = {0,0,0};
    borders.shader = shaders["curve"];

}



void scene_model::set_gui()
{
    // Can set the speed of the animation
    ImGui::SliderFloat("Time scale", &timer.scale, 0.05f, 2.0f, "%.2f s");
    ImGui::SliderFloat("Interval create sphere", &gui_scene.time_interval_new_sphere, 0.05f, 2.0f, "%.2f s");
    ImGui::Checkbox("Add sphere", &gui_scene.add_sphere);

    bool stop_anim  = ImGui::Button("Stop"); ImGui::SameLine();
    bool start_anim = ImGui::Button("Start");

    if(stop_anim)  timer.stop();
    if(start_anim) timer.start();
}





#endif
