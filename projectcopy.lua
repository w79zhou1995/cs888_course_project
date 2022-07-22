g1 = gr.material({0.13, 0.37, 0.105}, {0,0,0}, 0, 0)
g2 = gr.material({0.76, 0.8, 0.16}, {0,0,0}, 0, 0)
g3 = gr.material({0.38, 0.454, 0.02}, {0,0,0}, 0, 0)
g4 = gr.material({0.3, 0.7, 0.1}, {0,0,0}, 0, 0)
g5 = gr.material({0.51, 0.603, 0.305}, {0,0,0}, 0, 0)
tea = gr.material({0.47, 0.34, 0.21}, {0,0,0}, 20, 0.6) 
beige = gr.material({0.75, 0.73, 0.68}, {0,0,0}, 0, 0)
baby_blue = gr.material({0.72, 0.75, 0.81}, {0, 0, 0}, 0, 0)
jem1 = gr.material({1, 0.3, 0.0}, {0, 0, 0}, 10, 0.4)
butter_yellow = gr.material({0.94, 0.82, 0.61}, {0.0, 0.0, 0.0}, 2, 0)
bread_out = gr.material({0.73, 0.59, 0.51}, {0.0, 0.0, 0.0}, 0, 0)
bread_in = gr.material({0.94, 0.92, 0.84}, {0.0, 0.0, 0.0}, 0, 0)
china = gr.material({1, 1, 1}, {0.0, 0.0, 0.0}, 20, 0.2)
china_stripe = gr.material({-2, 1, 1}, {0.0, 0.0, 0.0}, 20, 0.2)
metal = gr.material({0.2, 0.2, 0.2}, {0.9, 0.9, 0.9}, 150, 0.8)
dark_blue = gr.material({0, 0, 0.4}, {0.0, 0.0, 0.0}, 0, 0)
blue = gr.material({0.4, 0.4, 1}, {0.0, 0.0, 0.0}, 0, 0)
black = gr.material({0, 0, 0}, {0.0, 0.0, 0.0}, 0, 0)
white = gr.material({1, 1, 1}, {0.0, 0.0, 0.0}, 0, 0)
pink = gr.material({1, 0.3, 0.3}, {0.0, 0.0, 0.0}, 0, 0)
skin = gr.material({1, 0.85, 0.85}, {0.0, 0.0, 0.0}, 0, 0)
stone = gr.material({0.8, 0.7, 0.7}, {0.0, 0.0, 0.0}, 0, 0)
grass = gr.material({0.1, 0.7, 0.1}, {0.0, 0.0, 0.0}, 0, 0)
mirror = gr.material({0.8, 0.7, 0.7}, {1.0, 1.0, 1.0}, 150, 1)
wood = gr.material({-1, 0, 0}, {0.0, 0.0, 0.0}, 0, 0)
wood_1 = gr.material({-3, 0, 0}, {0.0, 0.0, 0.0}, 0, 0)
cloth = gr.material({-4, 0, 0}, {0.0, 0.0, 0.0}, 0, 0)
grass = gr.material({0.1, 0.7, 0.1}, {0.0, 0.0, 0.0}, 0, 0)
hide = gr.material({0.84, 0.6, 0.53}, {0.3, 0.3, 0.3}, 20, 0.3)
-- ##############################################
-- the arc
-- ##############################################
-- ##############################################
-- the scene
-- ##############################################
time = 1005
start = 1000+5+5

--for i = 1, 10 do

scene = gr.node('scene')
scene:rotate('X', 23)

table_cloth1 = gr.cloth('tc1', 100, 'ball', {50, -20, 50}, 25, 0.1, time, start)
scene:add_child(table_cloth1)
table_cloth1:scale(0.064, 0.064, 0.064)
table_cloth1:translate(0, 5, 0)
table_cloth1:translate(-2.2, 4, 8.7)
table_cloth1:set_material(cloth)

s = gr.sphere('s')
table_cloth1:add_child(s)
s:set_material(hide)
s:scale(24.9, 24.9, 24.9)
s:translate(50, -20, 50)

name = 'test'.. tostring(time) .. '.png'

-- render it!
gr.render(scene,
	  name, 1000+5+5+5+5+5+5+5+5+5+5+5+5+5+5+5+5+5+5+5, 1000,
	  {0, 2, 30}, {0, 0, -1}, {0, 1, 0}, 50, {0.4, 0.4, 0.4}, 
{gr.light({200, 202, 430}, {0.8, 0.8, 0.8}, {1, 0, 0})})


--end
