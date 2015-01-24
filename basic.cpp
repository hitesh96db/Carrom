#include <iostream>
#include <cmath>
#include <stdlib.h>
#include <GL/glut.h>
#include <inttypes.h>
#include <stdio.h>
#include <sys/time.h>
#include <string.h>
#include <stdio.h>

#define PI 3.141592653589
#define WIDTH 0.05f

// ---------------------------------
// Game CONFIG
#define CONTROL 0.65f
#define RULE 0.4f
#define TIME_FACTOR 600
#define CALLS 20
#define VFACTOR 1.0f
#define SMASS 2
#define CMASS 1
#define TAKECARE 0.0f
#define NO_OF_COINS 7
#define SPACE 0.2f
#define FRICTION 0.2f
#define COLOR 1 // 1 = White, 0 = Black
#define WALL_COLLISION 1.5
#define ANGLE 0
 
// ---------------------------------
#define DEG2RAD(deg) (deg * PI / 180)
#define RAD2DEG(rad) (rad * 180 / PI)

char* itoa(int val, int base){
    
    static char buf[32] = {0};
    
    int i = 30;
    
    for(; val && i ; --i, val /= base)
    
        buf[i] = "0123456789abcdef"[val % base];
    
    return &buf[i+1];
    
}

using namespace std;
// ###########################################################################
// Coins

class Coin {
    public:
        Coin() {
            radius = 0.10f;
            velx = 0.0f;
            vely = 0.0f;
            stateX = 1;
            stateY = 1;
            rate = 0.0f;
            rest = 1;
            calculated = 0;
            last_check_time = -1;
            disappear = 0;
        }
        // Data Member
        char color;
        float posx,
              posy,
              velx,
              vely,
              radius,
              rate;
        double last_check_time;
        int stateX,
            stateY,
            rest,
            calculated,
            disappear;

        //Member Functions
        // Init
        void init(char c, float px, float py);
        void wallCollision();
        void coinCollision(int idx);
        void insidePocket();
};

void Coin::init(char c, float px, float py) {
    // Init
    color = c;
    posx = px;
    posy = py;
};

Coin c[7];

// ###########################################################################
// Pocket

class Pocket {
    public:
        // Data Member
        float posx;
        float posy;
        float radius;
        void init(float px, float py, float r);
};

void Pocket::init(float px, float py, float r) {
    posx = px;
    posy = py;
    radius = r;
};

Pocket pocket[4];

// ###########################################################################
// Striker

class Striker {
    // Variables
    public:
        Striker() {
            vely = 1.0f;
            velx = 1.0f;
            radius = 0.12f;
            powerStep = 1;
            stateX = 1;
            stateY = 1; 
            rate = 0.0f;
            found = 0;
            calculated = 0;
            last_check_time = -1;
        }
        float posx,
              posy,
              radius,
        // X-Y co-ordinates for helper line
              helperx,
              helpery,
              velx,
              vely,
              rate;
        // Power Step
        int powerStep,
            stateX,
            stateY,
            found,
            calculated;

        double last_check_time;

        // helper is the line for direction
        void init(float px, float py, float hx, float hy);
        void resetStates();
        void changeHelperRight();
        void changeHelperLeft();
        void increasePower();
        void decreasePower();
        void calculateVelocity();
        void wallCollision();
        void coinCollision();
        int insidePocket();
};

void Striker::init(float px, float py, float hx, float hy) {
    
    posx = px;
    posy = py;
    helperx = hx;
    helpery = hy;
};

Striker s;

// ###########################################################################
// Board : Outer Line + Inner Line

class Board {
    public:
        // Init
        Board() {
            olength = 5.0f;
            olineWidth = 20.0f; 
            ilength = 4.0f;
            ilineWidth = 3.0f;
            cradius = 0.55f;
            diff = 0.5f;
            pocketShift = 0.15f;
            pocketRadius = 0.15f;
            inAction = 0;
            gameOver = 0;
            // Circular Pattern Holes
        };

        // Data Member
        float olength,
              ilength,
              olineWidth,
              ilineWidth,
              cradius,
              diff,
              pocketShift,
              pocketRadius;
        int inAction, gameOver;

        // Member Functions
        void drawBoard();
        void drawPatterns();
        void drawCentreCircle();
        void drawCoins();
        void drawPockets();
        void drawStriker();
        void drawPower();
};

// Board Class Functions
void Board::drawBoard() {

    glColor3f(0.0f, 0.0f, 0.0f);
    
    // #############################
    // Draw Outer Board Line
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glLineWidth(olineWidth);
    glBegin(GL_QUADS);
    float half = olength/2;
    glVertex2f(-half, -half);
    glVertex2f(half, -half);
    glVertex2f(half, half);
    glVertex2f(-half, half);
    glEnd();
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    // #############################
    // Inner Line
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glColor3f(0.2f, 0.4f, 0.6f);
    // Draw Inner Line
    glLineWidth(ilineWidth);
    glBegin(GL_QUADS);
    half = ilength/2;
    glVertex2f(-half, -half);
    glVertex2f(half, -half);
    glVertex2f(half, half);
    glVertex2f(-half, half);
    glEnd();
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    glColor3f(0.2f, 0.4f, 0.6f);
    // Draw Inner Line
    glLineWidth(ilineWidth);
    glBegin(GL_QUADS);
    half = (ilength-diff)/2;
    glVertex2f(-half, -half);
    glVertex2f(half, -half);
    glVertex2f(half, half);
    glVertex2f(-half, half);
    glEnd();
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    
    int angle = 0;
    float dia = 2*cradius;
    glColor3f(0.2f, 0.4f, 0.6f);
    glPushMatrix();
    do {
        glLoadIdentity();
        glRotatef((float)angle,0.0,0.0,1.0);
        glTranslatef(0.0f,0.0f,-1.088f);
        glBegin(GL_LINE_LOOP);
        glVertex3f(0,0,0);
        glVertex3f(0,dia*0.10,0);
        glVertex3f(dia*0.01,dia*0.01,0);
        glEnd();
        // glColor3f(1.0f,1.0f,1.0f);
        glBegin(GL_POLYGON);
        glVertex3f(0,0,0);
        glVertex3f(0,dia*0.10,0);
        glVertex3f(dia*0.01,dia*0.01,0);
        glEnd();
        angle+=45;
    }while(angle!=360);
    glPopMatrix();
    glColor3f(0.0f, 0.0f, 0.0f);
    drawCentreCircle(); 
    drawPatterns();
    drawPockets();
    drawCoins();
    drawStriker();
    drawPower();
};

void Board::drawCentreCircle() {

    // ############################
    // Centre Circle    
    glBegin(GL_LINE_LOOP);
    float num_segments = 150;
    float theta, x, y, rad;
    rad = cradius;
    for (int ii = 0; ii < num_segments; ii++)   {
        theta = 2.0f * PI * float(ii) / float(num_segments);//get the current angle 
        x = rad * cosf(theta);//calculate the x component 
        y = rad * sinf(theta);//calculate the y component 
        glVertex2f(x, y);//output vertex 
    }
    glEnd();

};
void Board::drawPatterns() {

    // ############################
    // Pattern Holes

    float add_x, add_y, half;
    int i;
    float num_segments = 150;
    float theta, x, y, rad;

    half = ((ilength-(diff/2))/2);
    rad = diff/4;

    for(i=0;i<4;i++) {
        switch(i) {
            case 0:add_x = -half;
                   add_y = -half;
                   break;
            case 1:add_x = half;
                   add_y = -half;
                   break;
            case 2:add_x = half;
                   add_y = half;
                   break;
            case 3:add_x = -half;
                   add_y = half;
                   break;
        }
        glColor3f(1.0f, 0.0f, 0.0f);
        glBegin(GL_LINE_LOOP);
        for (int ii = 0; ii < num_segments; ii++)   {
            theta = 2.0f * PI * float(ii) / float(num_segments);//get the current angle 
            x = add_x + rad * cosf(theta);//calculate the x component 
            y = add_y + rad * sinf(theta);//calculate the y component 
            glVertex2f(x, y);//output vertex 
        }
        glEnd();
    }
};

void Board::drawPockets() {

    int j, i;
    float add_x, add_y, rad;

    for(i=0;i<4;i++) {

        rad = pocket[i].radius;
        add_x = pocket[i].posx;
        add_y = pocket[i].posy;
        glColor3f(0.0f, 0.0f, 0.0f);

        glBegin(GL_TRIANGLE_FAN);
        for(j=1; j<360 ; j++) {
            glVertex2f(add_x + rad * cos(DEG2RAD(j)), add_y + rad * sin(DEG2RAD(j)));
        }
        glEnd();
    }
};

void Board::drawCoins() {

    // ############################
    // Coins:
    char color;
    int j, i;
    float add_x, add_y, rad;

    for(i=0;i<NO_OF_COINS;i++) {

        if(c[i].disappear) 
            continue;

        rad = c[i].radius;
        add_x = c[i].posx;
        add_y = c[i].posy;
        color = c[i].color;
        if (color == 'r')
            glColor3f(1.0f, 0.0f, 0.0f);
        else if (color == 'b')
            glColor3f(0.0f, 0.0f, 0.0f);
        else
            glColor3f(1.0f, 1.0f, 1.0f);

        glBegin(GL_TRIANGLE_FAN);
        for(j=1; j<360 ; j++) {
            glVertex2f(add_x + rad * cos(DEG2RAD(j)), add_y + rad * sin(DEG2RAD(j)));
        }
        glEnd();
    }

};

void Board::drawStriker() {
    
    int j;
    float add_x, add_y, rad;
    rad = s.radius;
    add_x = s.posx;
    add_y = s.posy;
        
    glColor3f(0.58f, 0.098f, 0.098f);
    glBegin(GL_TRIANGLE_FAN);
    for(j=1; j<360 ; j++) {
        glVertex2f(add_x + rad * cos(DEG2RAD(j)), add_y + rad * sin(DEG2RAD(j)));
    }
    glEnd();

    if (!inAction) {
        glColor4f(0.2f, 0.4f, 0.6f, 0.25f);
        glBegin(GL_LINES);
        glVertex2f(add_x, add_y);
        glVertex2f(s.helperx, s.helpery);
        glEnd();
    }

};

void Board::drawPower() {
    
    float half = olength/2;
    float lx, ux, ly, uy;
    lx = ux = -half - 0.5f;
    ly = -half;
    uy = ly + s.powerStep*0.5f;

    glLineWidth(20.0f); 
    switch(s.powerStep) {
        case 5:glColor3f(1.0f, 0.0f, 0.0f);break;
        default:glColor3f(0.4, 0.8, 0.0f);
    } 
    
    glBegin(GL_LINES);
    glVertex2f(lx, ly);
    glVertex2f(ux, uy);
    glEnd();
};

Board board;

// ###########################################################################
// Game - stores physics stats

class Game {
    public:

        // Init
        Game() {
            force = FRICTION;
            scoreInt = 30;
        };
        double force;
        int scoreInt;
        int coinCol[7][7];
        int strikerCol[7];
        void init();
        void outputScore();
};

Game game;

void Game::init() {
int i, j;
    for(i=0;i<7;i++) {
        strikerCol[i] = 0;
        for(j=0;j<7;j++) {
            coinCol[i][j] = 0;
        }
    }
};

void Game::outputScore() {
   
    int sc = game.scoreInt;
    int minus = 0;
    if(sc < 0) {
        minus = 1;
        sc = -sc;
    }
    char *score = itoa(sc, 10);
    float half = board.olength/2;
    int len, i;
    len = (int)strlen(score);

    if( board.gameOver ) {
        glRasterPos3f(-board.cradius-0.3f, 0.3f, -8.0f);
        for(i=0;i<11;i++) {
            glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, "Game Over! "[i]);
        }
        glRasterPos3f(-board.cradius-0.3f, -0.3f, -8.0f);
        for(i=0;i<13;i++) {
            glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, "Final Score: "[i]);
        }
        glColor3f(1.0f, 0.0f, 0.0f);
        glRasterPos3f(-0.2f, -0.6f, -8.0f);
        if(minus) {
            glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, '-');
        }
        for(i=0;i<len;i++) {
            glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, score[i]);
        }
        return;
    }  
    glColor3f(0.0f, 0.0f, 0.0f);
    glRasterPos3f(-half-2.5f, half + 3.5f, -8.0f);
    for(i=0;i<7;i++) {
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, "Score: "[i]);
    }
    glColor3f(1.0f, 0.0f, 0.0f);
    glRasterPos3f(-half-1.0f, half + 3.5f, -8.0f);
    if(minus) {
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, '-');
    }
    for(i=0;i<len;i++) {
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, score[i]);
    }
}

void updateScore(int value) {
    if(!board.gameOver) {
        game.scoreInt--;
    }
    glutTimerFunc(1000, updateScore, 1);
};

// ###########################################################################
// Handler Functions

void drawScene(void) {

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glPushMatrix();
    glTranslatef(0.0f, 0.0f, -8.0f);
    board.drawBoard();
    // Display Score
    game.outputScore();

    glPopMatrix();
    glutSwapBuffers();

};
int dragflag = 1;

void checkGameOver() {

    char color = COLOR == 1 ? 'w': 'b';
    int i;
    board.gameOver = 1;
    for(i=0;i<7;i++) {
        if((c[i].color == color || c[i].color == 'r') && !c[i].disappear)
            board.gameOver = 0;
    }
}
int closeEnough(float spx, float spy, float sr, float cpx, float cpy, float cr, int takecare, int far=0) {

    float dist = sqrt((spx-cpx)*(spx-cpx) + (spy-cpy)*(spy-cpy));
    float add;
    if(takecare == 1) 
        add = TAKECARE;
    else
        add = 0.0f;
    if(far==1) {
        if(dist > sr + cr + add) 
            return 1;
        return 0;
    } 
    if( dist < sr + cr + add )
        return 1;
    return 0;
};


void fireCoin(int idx) {

   if(c[idx].disappear) 
        return;

    struct timeval tv;
    gettimeofday(&tv, NULL);
    double now = (tv.tv_sec) * 1000 + (tv.tv_usec) / 1000 ;
    double time_elapsed = 0;

    if (c[idx].last_check_time == -1) {
        c[idx].last_check_time = now;
    }
    else {
        time_elapsed = now - c[idx].last_check_time;
        c[idx].last_check_time = now;
    }
    // vx = ux + at
    // a
    float acc = game.force;
    // ux, uy
    float velx = c[idx].velx;
    float vely = c[idx].vely;
    time_elapsed /= TIME_FACTOR;
    // Y : v = u + at
    //   : s = ut + at2/2
    float new_vely = 0.0f;
    float disty = 0.0f;

    if(c[idx].stateY == 1)
        acc = -game.force;
    else
        acc = game.force;

    new_vely = vely + acc*time_elapsed;
    c[idx].vely = new_vely;
    vely /= VFACTOR;
    disty = vely*time_elapsed + ((acc*time_elapsed*time_elapsed)/2);

    if( c[idx].stateY == 1 ) {
        if(c[idx].vely <= 0.0f) {
            c[idx].vely = 0.0f;
            disty = 0.0f;
        }
    }
    else {
        if(c[idx].vely >= 0.0f) {
            c[idx].vely = 0.0f;
            disty = 0.0f;
        }
    }

    // --------------------------------------------------------------
    // X : v = u + at
    //   : s = ut + at2/2

    float new_velx = 0.0f;
    float distx = 0.0f;

    if( c[idx].stateX == 1 )
        acc = -game.force;
    else
        acc = game.force;

    new_velx = c[idx].rate*c[idx].vely;
    c[idx].velx = new_velx;
    velx /= VFACTOR;
    distx = velx*time_elapsed + ((acc*time_elapsed*time_elapsed)/2);

    if( c[idx].stateX == 1 ) {
        if(c[idx].velx <= 0.0f || c[idx].vely == 0.0f) {
            c[idx].velx = 0.0f;
            distx = 0.0f;
        }
    }
    else {
        if(c[idx].velx >= 0.0f || c[idx].vely == 0.0f) {
            c[idx].velx = 0.0f;
            distx = 0.0f;
        }
    }
    c[idx].posx += distx;
    c[idx].posy += disty;
    
    c[idx].wallCollision();
    c[idx].coinCollision(idx);
    c[idx].insidePocket();

    glutTimerFunc(CALLS, fireCoin, idx);

};

void coinCoinCollision(int idx1, int idx2) {

    // idx = coin index
    // Coin -------- 
    float cx1 = c[idx1].posx;
    float cy1 = c[idx1].posy;
    float cvx1 = c[idx1].velx;
    float cvy1 = c[idx1].vely;

    // Coin2 -------- 
    float cx2 = c[idx2].posx;
    float cy2 = c[idx2].posy;
    float cvx2 = c[idx2].velx;
    float cvy2 = c[idx2].vely;
    float angleLine, slope;
    //if(cx2-cx1 < 0) {
    //    angleLine = 0;
   // }
   // else {
    slope = ((float)(cy2 - cy1))/((float)(cx2 - cx1));
    angleLine = RAD2DEG(atanf(slope));
    if( angleLine < 0 ) {
        angleLine = 180 - abs(angleLine);
    }   
   // }

    float an = DEG2RAD(angleLine);
    float cos0 = cos(an);
    float sin0 = sin(an);

    float coin1_vx = cvx1*cos0 + cvy1*sin0;
    float coin1_vy = cvy1*cos0 - cvx1*sin0;

    float coin2_vx = cvx2*cos0 + cvy2*sin0;
    float coin2_vy = cvy2*cos0 - cvx2*sin0;

    // Collision
    float coin1_nvx = coin2_vx;
    float coin2_nvx = coin1_vx;

    // Coin1 New Velocity
    c[idx1].velx = (coin1_nvx*cos0 - coin1_vy*sin0);
    c[idx1].vely = (coin1_vy*cos0 + coin1_nvx*sin0);
    // Rate
    if(c[idx1].vely == 0.0f) 
        c[idx1].rate = 0.0f;
    else 
        c[idx1].rate = c[idx1].velx/c[idx1].vely;

    if(c[idx1].velx < 0.0f)
        c[idx1].stateX = -1;
    else
        c[idx1].stateX = 1;

    if(c[idx1].vely < 0.0f)
        c[idx1].stateY = -1;
    else
        c[idx1].stateY = 1;

    // Coin New Velocity
    c[idx2].velx = (coin2_nvx*cos0 - coin2_vy*sin0);
    c[idx2].vely = (coin2_vy*cos0 + coin2_nvx*sin0);
    // Rate
    if(c[idx2].vely == 0.0f)
        c[idx2].rate = 0.0f;
    else
        c[idx2].rate = c[idx2].velx/c[idx2].vely;

    if(c[idx2].velx < 0.0f)
        c[idx2].stateX = -1;
    else
        c[idx2].stateX = 1;

    if(c[idx2].vely < 0.0f)
        c[idx2].stateY = -1;
    else
        c[idx2].stateY = 1;
    
}

void strikerCoinCollision(int idx) {

    // idx = coin index
    // Coin -------- 
    float cx = c[idx].posx;
    float cy = c[idx].posy;
    float cvx = c[idx].velx;
    float cvy = c[idx].vely;
    
    // Striker -------- 
    float sx = s.posx;    
    float sy = s.posy;
    float svx = s.velx;
    float svy = s.vely;
    float slope, angleLine;
    // Line of Impact
    //if(sx-cx==0) {
    //    angleLine = 0;
   // }
    //else {
    slope = ((float)(sy - cy))/((float)(sx - cx));
    angleLine = RAD2DEG(atanf(slope));
    if( angleLine < 0 ) {
        angleLine = 180 - abs(angleLine);
    }
   // }
    
    float an = DEG2RAD(angleLine);
    float cos0 = cos(an);
    float sin0 = sin(an);
    
    float striker_vx = svx*cos0 + svy*sin0;
    float striker_vy = svy*cos0 - svx*sin0;
    
    float coin_vx = cvx*cos0 + cvy*sin0;
    float coin_vy = cvy*cos0 - cvx*sin0;
    
    // Collision
    float striker_nvx = ((striker_vx*(SMASS - CMASS) + 2*CMASS*coin_vx)/(CMASS + SMASS));
    float coin_nvx = ((coin_vx*(CMASS - SMASS) + 2*SMASS*striker_vx)/(CMASS + SMASS));
        
    // Striker New Velocity
    s.velx = (striker_nvx*cos0 - striker_vy*sin0);
    s.vely = (striker_vy*cos0 + striker_nvx*sin0);

    if(s.vely == 0.0f)
        s.rate = 0.0f;
    else 
        s.rate = s.velx/s.vely; 
    
    // Coin New Velocity
    c[idx].velx = (coin_nvx*cos0 - coin_vy*sin0);
    c[idx].vely = (coin_vy*cos0 + coin_nvx*sin0);

    if(c[idx].vely == 0.0f)
        c[idx].rate = 0.0f;
    else
        c[idx].rate = c[idx].velx/c[idx].vely;
    
    if(c[idx].velx < 0.0f) 
        c[idx].stateX = -1;
    else
        c[idx].stateX = 1;
    
    if(c[idx].vely < 0.0f) 
        c[idx].stateY = -1;
    else
        c[idx].stateY = 1;

    fireCoin(idx);
}

void Striker::resetStates() {
    vely = 1.0f;
    velx = 1.0f;
    radius = 0.12f;
    powerStep = 1;
    stateX = 1;
    stateY = 1;
    rate = 0.0f;
    found = 0;
    calculated = 0;
    last_check_time = -1;
}
void restartState() {
    
    board.inAction = 0;
    s.helperx = 0;
    s.helpery = (board.olength/2);
    // Set position of Striker
    float half = ((board.ilength - (board.diff/2))/2);
    s.init(0, -half, 0, (board.olength/2));
    s.resetStates();
    game.init();
};

void Coin::wallCollision() {

    float half = (board.olength/2);
    float rad = radius;
    if( (posx + rad + WIDTH >= half && stateX == 1) || (posx - rad - WIDTH <= -half && stateX == -1)) {
        stateX = -stateX;
        velx = -velx/WALL_COLLISION;
        vely = vely/WALL_COLLISION;
        rate = velx/vely;
    }
    if (( posy + rad + WIDTH >= half && stateY == 1)|| (posy - rad - WIDTH <= -half && stateY == -1 )) {
        stateY = -stateY;
        vely = -vely/WALL_COLLISION;
        velx = velx/WALL_COLLISION;
        rate = velx/vely;
    }
};

void Coin::insidePocket() {
    int i;
    if(disappear) 
        return;

    for(i=0;i<4;i++) {
        if(closeEnough(posx, posy, radius, pocket[i].posx, pocket[i].posy, pocket[i].radius, 0) == 1) {
            // Inside
            disappear = 1;
            char correct_color = COLOR == 1 ? 'w' : 'b';
            if(color == correct_color) 
                game.scoreInt += 10;
            else if(color == 'r')
                game.scoreInt += 50;
            else
                game.scoreInt -= 5;
            break;
        }
    }

};

void Coin::coinCollision(int idx) {

    if(disappear)
        return;

    int i;
    //float rate = velx/vely;
    int no_of_coins = NO_OF_COINS;
    for(i=0 ; i<no_of_coins ;i++) {
        // Check if collision takes place
        if(i!=idx && !c[i].disappear) {
            if((closeEnough(posx, posy, radius, c[i].posx, c[i].posy, c[i].radius, 1) == 1) && (!game.coinCol[idx][i])) {
                // Collide !
                // Pass coin Index
                coinCoinCollision(idx, i);
                game.coinCol[idx][i] = 1;
                game.coinCol[i][idx] = 1;
            }
            else if((closeEnough(posx, posy, radius, c[i].posx, c[i].posy, c[i].radius, 1, 1) == 1) && game.coinCol[idx][i]) {
                game.coinCol[idx][i] = 0;
                game.coinCol[i][idx] = 0;
            }
        }
    }
};

void Striker::wallCollision() {

    float half = (board.olength/2);
    float rad = radius;
    if( (posx + rad + WIDTH >= half && stateX == 1) || (posx - rad - WIDTH <= -half && stateX == -1)) {
        stateX = -stateX;
        velx = -velx/WALL_COLLISION;
        vely = vely/WALL_COLLISION;
        rate = velx/vely;
    }
    if (( posy + rad + WIDTH >= half && stateY == 1)|| (posy - rad - WIDTH <= -half && stateY == -1 )) {
        stateY = -stateY;
        vely = -vely/WALL_COLLISION;
        velx = velx/WALL_COLLISION;
        rate = velx/vely;
    }
};

int Striker::insidePocket() {
    int i;
    for(i=0;i<4;i++) {
        if(closeEnough(posx, posy, radius, pocket[i].posx, pocket[i].posy, pocket[i].radius, 0) == 1) {
            // Inside
            game.scoreInt -= 5;
            restartState();
            return 1;
        }
    }
    return 0;
};

void Striker::coinCollision() {

    int i;
    int no_of_coins = NO_OF_COINS;
    for(i=0 ; i<no_of_coins ;i++) {
        // Check if collision takes place
        if(c[i].disappear)
            continue; 
        if((closeEnough(posx, posy, radius, c[i].posx, c[i].posy, c[i].radius, 1) == 1)&&(!game.strikerCol[i])) {
            // Collide !
            // Pass coin Index
            strikerCoinCollision(i);
            game.strikerCol[i] = 1;
        }
        else if((closeEnough(posx, posy, radius, c[i].posx, c[i].posy, c[i].radius, 1, 1) == 1)&&(game.strikerCol[i]))
            game.strikerCol[i] = 0;
    }
};

void Striker::calculateVelocity() {

    calculated = 1;
    float x1 = posx;
    float y1 = posy;
    float x2 = helperx;
    float y2 = helpery;
    float slope = (y2 - y1)/(x2 - x1);
    float angle = RAD2DEG(atan(slope));
    if ( angle < 0.0f ) {
        stateX = -1;
    }
    else {
        stateX = 1;
    }
    float sin0 = sin(DEG2RAD(abs(angle)));
    float cos0 = cos(DEG2RAD(abs(angle)));
    vely = powerStep * sin0 * CONTROL;
    velx = stateX * powerStep * cos0 * CONTROL;
    rate = s.velx/s.vely;
}

void Striker::changeHelperLeft() {
    // Move the helper line to the left
    float hx = helperx;
    float hy = helpery;
    float half = (board.olength/2);
    float half2 = ((board.ilength-board.diff)/2);

    if( hy == half ) {
        helperx -= 0.1f;
        if (helperx < -half) {
            helperx = -half;
            helpery = half - 0.1f;
        }
    }
    if( hx == -half ) {
        helpery -= 0.1f;
        if (helpery < -half2 + RULE) {
            helpery = -half2 + RULE;
        }
    }
    if( hx == half ) {
        helpery += 0.1f;
        if (helpery > half) {
            helperx = half-0.1f;
            helpery = half;
        }
    }

    calculateVelocity();
};

void Striker::changeHelperRight() {
    float hx = helperx;
    float hy = helpery;
    float half = (board.olength/2);
    float half2 = ((board.ilength-board.diff)/2);

    if( hy == half ) {
        helperx += 0.1f;
        if (helperx > half) {
            helperx = half;
            helpery = half - 0.1f;
        }
    }
    if( hx == -half ) {
        helpery += 0.1f;
        if (helpery > half) {
            helpery = half;
            helperx = -half+0.1f;
        }
    }
    if( hx == half ) {
        helpery -= 0.1f;
        if (helpery < -half2 + RULE) {
            helpery = -half2 + RULE;
        }
    } 
    calculateVelocity();
};

void Striker::increasePower() {
    powerStep += 1;
    calculateVelocity();
    if (powerStep >= 5) {
        powerStep = 5;
    }
};  

void Striker::decreasePower() {
    powerStep -= 1;
    calculateVelocity();
    if(powerStep < 0 ) {
        powerStep = 0;
    }
};

void fireStriker(int value) {
    
    if(!board.inAction) 
        return;

    if(!s.calculated) {
        // Find Initial Velocity if not found
        s.calculateVelocity();
    }

    struct timeval tv;
    gettimeofday(&tv, NULL);
    double now = (tv.tv_sec) * 1000 + (tv.tv_usec) / 1000 ;
    double time_elapsed = 0;

    if (s.last_check_time == -1) {
        s.last_check_time = now;
    }
    else {
        time_elapsed = now - s.last_check_time;
        s.last_check_time = now;
    }   
    // vx = ux + at
    // a
    float acc = game.force;
    // ux, uy
    float velx = s.velx;
    float vely = s.vely;
    time_elapsed /= TIME_FACTOR;
    
    // Y : v = u + at
    //   : s = ut + at2/2
    float new_vely = 0.0f;
    float disty = 0.0f;

    if(s.stateY == 1)
        acc = -game.force;
    else 
        acc = game.force;
    
    new_vely = vely + acc*time_elapsed;
    s.vely = new_vely;
    vely /= VFACTOR;
    disty = vely*time_elapsed + ((acc*time_elapsed*time_elapsed)/2);
    if( s.stateY == 1 ) {
        if(s.vely <= 0.0f) {
            s.vely = 0.0f;
            disty = 0.0f;
        }
    }   
    else {
        if(s.vely >= 0.0f) {
            s.vely = 0.0f;
            disty = 0.0f;
        }
    }
 
    // --------------------------------------------------------------
    // X : v = u + at
    //   : s = ut + at2/2

    float new_velx = 0.0f;
    float distx = 0.0f;

    if( s.stateX == 1 )
        acc = -game.force;
    else
        acc = game.force;

    new_velx = s.rate*s.vely;
    s.velx = new_velx;
    velx /= VFACTOR;
    distx = velx*time_elapsed + ((acc*time_elapsed*time_elapsed)/2);

    if( s.stateX == 1 ) {
        if(s.velx <= 0.0f) {
            s.velx = 0.0f;
            distx = 0.0f;
        }
    }
    else {
        if(s.velx >= 0.0f) {
            s.velx = 0.0f;
            distx = 0.0f;
        }
    }
    
    s.posx += distx;
    s.posy += disty;

    s.wallCollision();
    s.coinCollision();
    if(s.insidePocket()) {
        return;
    }
    checkGameOver();
    glutTimerFunc(CALLS, fireStriker, 1);

}


void handleClick(int px, int py) {

    int w = glutGet(GLUT_SCREEN_WIDTH);
    float x = px - w / 2.0;
    float half = board.ilength/2;
    if(dragflag && (abs(x) < 200)) {
        s.posx = (x / 200) * half * 0.95;
        half = ((board.ilength - (board.diff/2))/2);
        s.posy = -half;
    }
};

// Called when click is pressed
void setFlag(int button, int state, int x, int y) {

    if(button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN)
        dragflag = 1;
    else
        dragflag = 0;

    if(button == GLUT_LEFT_BUTTON) {
        if(!board.inAction) {
            board.inAction = 1;
            glutTimerFunc(10, fireStriker, 1);
        }
    }
};


void handleHover(int px, int py) {

    int w = glutGet(GLUT_SCREEN_WIDTH);
    int h = glutGet(GLUT_SCREEN_HEIGHT);
    float x = px - w / 2.0;
    float y = (float)(((float)h / 2.0) - (float)py);
    float half = board.ilength/2;
    s.helperx = (x / 200) * half * 0.8;
    s.helpery = (y / 200) * half * 0.8;
    float dist = sqrt((s.posx-s.helperx)*(s.posx-s.helperx) + (s.posy-s.helpery)*(s.posy-s.helpery));
    if(!board.inAction) {
        s.powerStep = ceil(dist);
        if(s.powerStep > 5)
            s.powerStep = 5;
    }
    if(s.helpery < s.posy ) {
        s.helpery = s.posy + 0.5f;
    }
    glutPostRedisplay();
};


void handleArrowKeys(int key, int x, int y) {

    if(board.inAction)
        return;

    float half2 = (board.ilength)/2;
    if (key == GLUT_KEY_LEFT) {
        s.posx -= 0.1f;
        if ( s.posx - s.radius < -half2 ) {
            s.posx = -half2 + s.radius;
        }
    }
    else if(key == GLUT_KEY_RIGHT) {
        s.posx += 0.1f;
        if ( s.posx + s.radius > half2 ) {
            s.posx = half2 - s.radius;
        }
    }
    else if(key == GLUT_KEY_UP) {
        s.increasePower();
    }
    else if(key == GLUT_KEY_DOWN) {
        s.decreasePower();
    }
};

void GameSetup() {

    game.init();
    // Put all coins into the centre

    // Red
    char color;
    float posx,
          posy;

    float rad = 2*c[0].radius + SPACE;

    color = 'r';
    posx = 0.0f;
    posy = 0.0f;
    c[0].init(color, posx, posy);

    // White 1 - 3
    color = 'w';
    float deg = 270 + ANGLE;

    posx = rad * cos(DEG2RAD(deg));
    posy = rad * sin(DEG2RAD(deg));
    c[1].init(color, posx, posy);

    // + 120
    deg = deg + 120;
    posx = rad * cos(DEG2RAD(deg));
    posy = rad * sin(DEG2RAD(deg));
    c[2].init(color, posx, posy);

    // + 120
    deg = deg + 120;
    posx = rad * cos(DEG2RAD(deg));
    posy = rad * sin(DEG2RAD(deg));
    c[3].init(color, posx, posy);

    // Black 4 - 6
    color = 'b';

    deg = 330 + ANGLE;
    posx = rad * cos(DEG2RAD(deg));
    posy = rad * sin(DEG2RAD(deg));
    c[4].init(color, posx, posy);

    deg = deg + 120;
    posx = rad * cos(DEG2RAD(deg));
    posy = rad * sin(DEG2RAD(deg));
    c[5].init(color, posx, posy);

    deg = deg + 120;
    posx = rad * cos(DEG2RAD(deg));
    posy = rad * sin(DEG2RAD(deg));
    c[6].init(color, posx, posy);

    // Set position of Pockets
    float half = board.olength/2;
    float shift = board.pocketShift;
    float prad = board.pocketRadius;
    pocket[0].init(-half+shift, -half+shift, prad);
    pocket[1].init(half-shift, -half+shift, prad);
    pocket[2].init(half-shift, half-shift, prad);
    pocket[3].init(-half+shift, half-shift, prad);

    // Set position of Striker
    half = ((board.ilength - (board.diff/2))/2);
    s.init(0, -half, 0, (board.olength/2));
    int i;
    for(i=0;i<NO_OF_COINS;i++) {
        glutTimerFunc(CALLS, fireCoin, i);
    }
    glutTimerFunc(1000, updateScore, 1);
};

void handleKeyPress(unsigned char key, int x, int y) {

    if((board.inAction && key!=114) || board.gameOver) {
        if( key != 27 )
            return;
    }

    switch(key) {
        case 27: exit(0);
                 break;
        case 97: s.changeHelperLeft();
                 break;
        case 99: s.changeHelperRight();
                 break;
        case 32: board.inAction = 1;
                 fireStriker(1);
                 break;
        case 114: board.inAction = 0; 
                  restartState();
                  break;
    }
};

void handleResize(int w, int h) {

    // Handle Resize
    glViewport(0, 0, w, h);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(45.0f, (float)w / (float)h, 1.0f, 200.0f);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

};

void GLSetup() {

    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
    int w = glutGet(GLUT_SCREEN_WIDTH);
    int h = glutGet(GLUT_SCREEN_HEIGHT);
    //int windowWidth = w * 3 / 4;
    //int windowHeight = h * 3 / 4;

    glutInitWindowSize(w, h);
    //glutInitWindowPosition((w - windowWidth) / 2, (h - windowHeight) / 2);
    glutCreateWindow("- Carrom -");
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    // Board Color
    glClearColor(1.0f, 0.82f, 0.6f, 0.0f);
    // Set Handler Functions
    glutDisplayFunc(drawScene);
    glutKeyboardFunc(handleKeyPress);
    glutSpecialFunc(handleArrowKeys);
    glutIdleFunc(drawScene);
    glutReshapeFunc(handleResize);
    glutPassiveMotionFunc(handleHover);
    glutMouseFunc(setFlag);
    glutMotionFunc(handleClick);
};

int main(int argc, char **argv) {
    
    glutInit(&argc, argv);
    GLSetup();
    GameSetup();
    // Main loop
    glutMainLoop();

return 0;
}
