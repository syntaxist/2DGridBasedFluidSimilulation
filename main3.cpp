#include <iostream>
#include <SFML/Graphics.hpp>



int resolution = 800;

int numberOfGridRow = 128;
int N = numberOfGridRow ;
float cellSize = (float)resolution / (float)numberOfGridRow;
const int n = 128 * 128;
int size = n;
int iter = 16;
float dt = 1;
float diff = 0.0001;
float visc = 0.0001;

float s[n];
float density[n];

float Vx[n];
float Vy[n];

float Vx0[n];
float Vy0[n];




int IX(int x, int y) {
    if (x>=0 &&  x<= N - 1 && y >= 0 && y <= N - 1)
    {
        return x + (y * N);
    }

}

void addDensity(int x, int y, float amount) {

    int index = IX((int)(std::floor(x / cellSize)), (int)(std::floor(y / cellSize)));
    density[index] += amount;
}

void addVelocity(int x, int y, float amountX, float amountY) {
    int index = IX((int)(std::floor(x / cellSize)), (int)(std::floor(y / cellSize)));
    Vx[index] += amountX;
    Vy[index] += amountY;
}


void set_bnd(int b, float x[]) {
    for (int i = 1; i < N - 1; i++) {
        x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
        x[IX(i, N - 1)] = b == 2 ? -x[IX(i, N - 2)] : x[IX(i, N - 2)];
    }
    for (int j = 1; j < N - 1; j++) {
        x[IX(0, j)] = b == 1 ? -x[IX(1, j)] : x[IX(1, j)];
        x[IX(N - 1, j)] = b == 1 ? -x[IX(N - 2, j)] : x[IX(N - 2, j)];
    }

    x[IX(0, 0)] = 0.5f * (x[IX(1, 0)] + x[IX(0, 1)]);
    x[IX(0, N - 1)] = 0.5f * (x[IX(1, N - 1)] + x[IX(0, N - 2)]);
    x[IX(N - 1, 0)] = 0.5f * (x[IX(N - 2, 0)] + x[IX(N - 1, 1)]);
    x[IX(N - 1, N - 1)] = 0.5f * (x[IX(N - 2, N - 1)] + x[IX(N - 1, N - 2)]);
}

void lin_solve(int b, float x[], float x0[], float a, float c) {
    float cRecip = 1.0 / c;
    for (int k = 0; k < iter; k++) {
        for (int j = 1; j < N - 1; j++) {
            for (int i = 1; i < N - 1; i++) {
                x[IX(i, j)] = (x0[IX(i, j)] + a * (x[IX(i + 1, j)] + x[IX(i - 1, j)] + x[IX(i, j + 1)] + x[IX(i, j - 1)])) * cRecip;
            }
        }
        set_bnd(b, x);
    }
}


void diffuse(int b, float x[], float x0[], float diff, float dt) {
    float a = dt * diff * (N - 2) * (N - 2);
    lin_solve(b, x, x0, a, 1 + 4 * a);
}


void project(float velocX[], float velocY[], float p[], float div[]) {
    for (int j = 1; j < N - 1; j++) {
        for (int i = 1; i < N - 1; i++) {
            div[IX(i, j)] = -0.5f * (
                velocX[IX(i + 1, j)]
                - velocX[IX(i - 1, j)]
                + velocY[IX(i, j + 1)]
                - velocY[IX(i, j - 1)]
                ) / N;
            p[IX(i, j)] = 0;
        }
    }

    set_bnd(0, div);
    set_bnd(0, p);
    lin_solve(0, p, div, 1, 4);

    for (int j = 1; j < N - 1; j++) {
        for (int i = 1; i < N - 1; i++) {
            velocX[IX(i, j)] -= 0.5f * (p[IX(i + 1, j)]
                - p[IX(i - 1, j)]) * N;
            velocY[IX(i, j)] -= 0.5f * (p[IX(i, j + 1)]
                - p[IX(i, j - 1)]) * N;
        }
    }
    set_bnd(1, velocX);
    set_bnd(2, velocY);
}


void advect(int b, float d[], float d0[], float velocX[], float velocY[], float dt) {
    float i0, i1, j0, j1;

    float dtx = dt * (N - 2);
    float dty = dt * (N - 2);

    float s0, s1, t0, t1;
    float tmp1, tmp2, x, y;

    float Nfloat = N;
    float ifloat, jfloat;
    int i, j;

    for (j = 1, jfloat = 1; j < N - 1; j++, jfloat++) {
        for (i = 1, ifloat = 1; i < N - 1; i++, ifloat++) {
            tmp1 = dtx * velocX[IX(i, j)];
            tmp2 = dty * velocY[IX(i, j)];
            x = ifloat - tmp1;
            y = jfloat - tmp2;

            if (x < 0.5f) x = 0.5f;
            if (x > Nfloat + 0.5f) x = Nfloat + 0.5f;
            i0 = std::floor(x);
            i1 = i0 + 1.0f;
            if (y < 0.5f) y = 0.5f;
            if (y > Nfloat + 0.5f) y = Nfloat + 0.5f;
            j0 = std::floor(y);
            j1 = j0 + 1.0f;

            s1 = x - i0;
            s0 = 1.0f - s1;
            t1 = y - j0;
            t0 = 1.0f - t1;

            int i0i = int(i0);
            int i1i = int(i1);
            int j0i = int(j0);
            int j1i = int(j1);

            d[IX(i, j)] =
                s0 * (t0 * d0[IX(i0i, j0i)] + t1 * d0[IX(i0i, j1i)]) +
                s1 * (t0 * d0[IX(i1i, j0i)] + t1 * d0[IX(i1i, j1i)]);
        }
    }

    set_bnd(b, d);
}

void step() {

    diffuse(1, Vx0, Vx, visc, dt);
    diffuse(2, Vy0, Vy, visc, dt);

    project(Vx0, Vy0, Vx, Vy);

    advect(1, Vx, Vx0, Vx0, Vy0, dt);
    advect(2, Vy, Vy0, Vx0, Vy0, dt);

    project(Vx, Vy, Vx0, Vy0);

    diffuse(0, s, density, diff, dt);
    advect(0, density, s, Vx, Vy, dt);
}


sf::VertexArray newShape(int numberOfQuads) {
    sf::VertexArray shape1(sf::Quads);
    for (size_t j = 0; j < numberOfQuads; j++)
    {
        for (size_t i = 0; i < numberOfQuads; i++)
        {
            sf::Vertex v1(sf::Vector2f(i * cellSize, j * cellSize));
            sf::Vertex v2(sf::Vector2f(v1.position.x, v1.position.y + cellSize));
            sf::Vertex v3(sf::Vector2f(v1.position.x + cellSize, v1.position.y));
            sf::Vertex v4(sf::Vector2f(v1.position.x + cellSize, v1.position.y + cellSize));

            v1.color = v2.color = v3.color = v4.color = sf::Color(255, 255, 255, rand() % 255);

            shape1.append(v2);
            shape1.append(v4);
            shape1.append(v3);
            shape1.append(v1);
        }
    }
    return shape1;
}

void changeColor(sf::VertexArray& array, float theDens[]) {
    for (int i = 0; i < N + 2; i++)
    {
        for (int j = 0; j < N + 2; j++)
        {
            int index = IX(i, j);
            index = index * 4;
            array[index].color = 
            array[index + 1].color = 
            array[index + 2].color = 
            array[index + 3].color = sf::Color(255, 255, 255, theDens[IX(i, j)] *100);
        }

    }

}


int main() {

    sf::Clock dtClock, fpsTimer;
    sf::RenderWindow window(sf::VideoMode(resolution, resolution), "Too Slow");
    sf::VertexArray shapes = newShape(numberOfGridRow);

    sf::Vector2f mousePressed;
    sf::Vector2f mouseReleased;
    bool mouseLeftPressed = false;
    window.setFramerateLimit(240);

    sf::Vector2f pMouse;
    sf::Vector2f cMouse;
    while (window.isOpen()) {




        float dtFR = dtClock.restart().asSeconds();
        if (fpsTimer.getElapsedTime().asSeconds() > 1) {
            fpsTimer.restart();
        }
        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                window.close();
            if (event.type == sf::Event::MouseButtonPressed)
            {
                if (event.mouseButton.button == sf::Mouse::Left)
                {
                    mouseLeftPressed = true;

                }
            }
            if (event.type == sf::Event::MouseButtonReleased)
            {
                if (event.mouseButton.button == sf::Mouse::Left)
                {
                    mouseLeftPressed = false;

                }
            }
            if (mouseLeftPressed == true)
            {
                if (event.type == sf::Event::MouseMoved)
                {

                    pMouse.x = sf::Mouse::getPosition(window).x;
                    pMouse.y = sf::Mouse::getPosition(window).y;
                    sf::sleep(sf::seconds(dtFR));
                    cMouse.x = sf::Mouse::getPosition(window).x;
                    cMouse.y = sf::Mouse::getPosition(window).y;
                    std::cout << "new mouse x: " << (cMouse.x - pMouse.x) << "," << (cMouse.y - pMouse.y) << std::endl;

                    addDensity(cMouse.x, cMouse.y, 100);
                    addVelocity(cMouse.x, cMouse.y, cMouse.x - pMouse.x, cMouse.y - pMouse.y);

                    window.clear(sf::Color(0, 0, 0));

                    step();
                    changeColor(shapes, density);


                    window.draw(shapes);
                    window.display();


                }
            }

        }
        window.clear(sf::Color(0, 0, 0));

        step();
        changeColor(shapes, density);


        window.draw(shapes);
        window.display();


    }
}
