%CMPEN 454 Project 2

%Load the necessary data files
load('vue2CalibInfo.mat')
load('vue4CalibInfo.mat')
%Initialization for Vue2
filenamevue2mp4 = 'Subject4-Session3-24form-Full-Take4-Vue2.mp4';
vue2video = VideoReader(filenamevue2mp4);
%Initialization for Vue4
filenamevue4mp4 = 'Subject4-Session3-24form-Full-Take4-Vue4.mp4';
vue4video = VideoReader(filenamevue2mp4);

%Get the world coordinates(x,y,z) and confidence values for each joint in
%every frame
frame = 1:1:26214;
right_shoulder_world = get_joint_world_coords(frame,1);
right_elbow_world = get_joint_world_coords(frame,2);
right_wrist_world = get_joint_world_coords(frame,3);
left_shoulder_world = get_joint_world_coords(frame,4);
left_elbow_world = get_joint_world_coords(frame,5);
left_wrist_world = get_joint_world_coords(frame,6);
right_hip_world = get_joint_world_coords(frame,7);
right_knee_world = get_joint_world_coords(frame,8);
right_ankle_world = get_joint_world_coords(frame,9);
left_hip_world = get_joint_world_coords(frame,10);
left_knee_world = get_joint_world_coords(frame,11);
left_ankle_world = get_joint_world_coords(frame,12);

%vue2 camera parameters
vue2_Rmat = [-0.759251863679378,-0.649110013194741,0.0468273239649182,0;-0.138093351133920,0.0903736400713044,-0.986287398050318,0;0.635977070240741,-0.755307087184137,-0.158254131628272,0; 0,0,0,1];
vue2_Tmat = [1,0,0,4450.1; 0,1,0,-5557.9; 0,0,1,-1949.1;0,0,0,1];
%vue2_transposed_rmat = transpose(vue2.Rmat);
%vue2_t_vec = [137.715365606944;805.522696238255;7336.52052400539];

%vue4 camera parameters
vue4_Rmat = [-0.816427158217721,0.576669940506233,0.0299745732425930,0;0.103963925393142,0.197851965453977,-0.974703084012177,0;-0.568012497698235,-0.792657794689825,-0.221484588574639,0; 0,0,0,1];
vue4_Tmat = [1,0,0,-4423.56297868868; 0,1,0,-5490.36044079971; 0,0,1,-1889.01674274890;0,0,0,1];
%vue4_transposed_rmat = transpose(vue4.Rmat);
%vue4_t_vec = [388.768652422159;295.060869156745;7283.00415138530];

%Get the pixel values for vue4 in vue2 and for vue2 in vue4
vue2_coords_in_vue4 = get_v2_in_v4();
vue4_coords_in_vue2 = get_v4_in_v2();
%Get the pixel values for all joints in Vue2
right_shoulder_pixel_vue2 = get_vue2_joint_pixel_coords(frame,right_shoulder_world);
right_elbow_pixel_vue2 = get_vue2_joint_pixel_coords(frame,right_elbow_world);
right_wrist_pixel_vue2 = get_vue2_joint_pixel_coords(frame,right_wrist_world);
left_shoulder_pixel_vue2 = get_vue2_joint_pixel_coords(frame,left_shoulder_world);
left_elbow_pixel_vue2 = get_vue2_joint_pixel_coords(frame,left_elbow_world);
left_wrist_pixel_vue2 = get_vue2_joint_pixel_coords(frame,left_wrist_world);
right_hip_pixel_vue2 = get_vue2_joint_pixel_coords(frame,right_hip_world);
right_knee_pixel_vue2 = get_vue2_joint_pixel_coords(frame,right_knee_world);
right_ankle_pixel_vue2 = get_vue2_joint_pixel_coords(frame,right_ankle_world);
left_hip_pixel_vue2 = get_vue2_joint_pixel_coords(frame,left_hip_world);
left_knee_pixel_vue2 = get_vue2_joint_pixel_coords(frame,left_knee_world);
left_ankle_pixel_vue2 = get_vue2_joint_pixel_coords(frame,left_ankle_world);
%Get the pixel values for all joints in Vue4
right_shoulder_pixel_vue4 = get_vue4_joint_pixel_coords(frame,right_shoulder_world);
right_elbow_pixel_vue4 = get_vue4_joint_pixel_coords(frame,right_elbow_world);
right_wrist_pixel_vue4 = get_vue4_joint_pixel_coords(frame,right_wrist_world);
left_shoulder_pixel_vue4 = get_vue4_joint_pixel_coords(frame,left_shoulder_world);
left_elbow_pixel_vue4 = get_vue4_joint_pixel_coords(frame,left_elbow_world);
left_wrist_pixel_vue4 = get_vue4_joint_pixel_coords(frame,left_wrist_world);
right_hip_pixel_vue4 = get_vue4_joint_pixel_coords(frame,right_hip_world);
right_knee_pixel_vue4 = get_vue4_joint_pixel_coords(frame,right_knee_world);
right_ankle_pixel_vue4 = get_vue4_joint_pixel_coords(frame,right_ankle_world);
left_hip_pixel_vue4 = get_vue4_joint_pixel_coords(frame,left_hip_world);
left_knee_pixel_vue4 = get_vue4_joint_pixel_coords(frame,left_knee_world);
left_ankle_pixel_vue4 = get_vue4_joint_pixel_coords(frame,left_ankle_world);

%4. Reconstruct the 3D location of each joint in the world coordinate system from the
%projected 2D joints you produced in Step3, using two-camera triangulation.

test_rec = calc_triangulation(1, left_wrist_pixel_vue4);

%5. Compute Euclidean (L²) distance between all joint pairs. This is a per joint, per frame L²
%distance between the original 3D joints and the reconstructed 3D joints providing a
%quantitative analysis of the distance between the joint pairs.



%Plotting to verify results:
%The following code creates images for selected frames and plots the joints
%and a created skeleton for vue4
figure(1)
verify_3Dto2D_vue2(1000, right_shoulder_pixel_vue2, right_elbow_pixel_vue2, right_wrist_pixel_vue2, left_shoulder_pixel_vue2, left_elbow_pixel_vue2, left_wrist_pixel_vue2, right_hip_pixel_vue2, right_knee_pixel_vue2, right_ankle_pixel_vue2, left_hip_pixel_vue2, left_knee_pixel_vue2, left_ankle_pixel_vue2)
figure(2)
verify_3Dto2D_vue2(5000, right_shoulder_pixel_vue2, right_elbow_pixel_vue2, right_wrist_pixel_vue2, left_shoulder_pixel_vue2, left_elbow_pixel_vue2, left_wrist_pixel_vue2, right_hip_pixel_vue2, right_knee_pixel_vue2, right_ankle_pixel_vue2, left_hip_pixel_vue2, left_knee_pixel_vue2, left_ankle_pixel_vue2)
figure(3)
verify_3Dto2D_vue2(10000, right_shoulder_pixel_vue2, right_elbow_pixel_vue2, right_wrist_pixel_vue2, left_shoulder_pixel_vue2, left_elbow_pixel_vue2, left_wrist_pixel_vue2, right_hip_pixel_vue2, right_knee_pixel_vue2, right_ankle_pixel_vue2, left_hip_pixel_vue2, left_knee_pixel_vue2, left_ankle_pixel_vue2)
figure(4)
verify_3Dto2D_vue2(15000, right_shoulder_pixel_vue2, right_elbow_pixel_vue2, right_wrist_pixel_vue2, left_shoulder_pixel_vue2, left_elbow_pixel_vue2, left_wrist_pixel_vue2, right_hip_pixel_vue2, right_knee_pixel_vue2, right_ankle_pixel_vue2, left_hip_pixel_vue2, left_knee_pixel_vue2, left_ankle_pixel_vue2)
figure(5)
verify_3Dto2D_vue2(20000, right_shoulder_pixel_vue2, right_elbow_pixel_vue2, right_wrist_pixel_vue2, left_shoulder_pixel_vue2, left_elbow_pixel_vue2, left_wrist_pixel_vue2, right_hip_pixel_vue2, right_knee_pixel_vue2, right_ankle_pixel_vue2, left_hip_pixel_vue2, left_knee_pixel_vue2, left_ankle_pixel_vue2)
%The following code creates images for selected frames and plots the joints
%and a created skeleton for vue4
figure(6)
verify_3Dto2D_vue4(1000, right_shoulder_pixel_vue4, right_elbow_pixel_vue4, right_wrist_pixel_vue4, left_shoulder_pixel_vue4, left_elbow_pixel_vue4, left_wrist_pixel_vue4, right_hip_pixel_vue4, right_knee_pixel_vue4, right_ankle_pixel_vue4, left_hip_pixel_vue4, left_knee_pixel_vue4, left_ankle_pixel_vue4)
figure(7)
verify_3Dto2D_vue4(5000, right_shoulder_pixel_vue4, right_elbow_pixel_vue4, right_wrist_pixel_vue4, left_shoulder_pixel_vue4, left_elbow_pixel_vue4, left_wrist_pixel_vue4, right_hip_pixel_vue4, right_knee_pixel_vue4, right_ankle_pixel_vue4, left_hip_pixel_vue4, left_knee_pixel_vue4, left_ankle_pixel_vue4)
figure(8)
verify_3Dto2D_vue4(10000, right_shoulder_pixel_vue4, right_elbow_pixel_vue4, right_wrist_pixel_vue4, left_shoulder_pixel_vue4, left_elbow_pixel_vue4, left_wrist_pixel_vue4, right_hip_pixel_vue4, right_knee_pixel_vue4, right_ankle_pixel_vue4, left_hip_pixel_vue4, left_knee_pixel_vue4, left_ankle_pixel_vue4)
figure(9)
verify_3Dto2D_vue4(15000, right_shoulder_pixel_vue4, right_elbow_pixel_vue4, right_wrist_pixel_vue4, left_shoulder_pixel_vue4, left_elbow_pixel_vue4, left_wrist_pixel_vue4, right_hip_pixel_vue4, right_knee_pixel_vue4, right_ankle_pixel_vue4, left_hip_pixel_vue4, left_knee_pixel_vue4, left_ankle_pixel_vue4)
figure(10)
verify_3Dto2D_vue4(20000, right_shoulder_pixel_vue4, right_elbow_pixel_vue4, right_wrist_pixel_vue4, left_shoulder_pixel_vue4, left_elbow_pixel_vue4, left_wrist_pixel_vue4, right_hip_pixel_vue4, right_knee_pixel_vue4, right_ankle_pixel_vue4, left_hip_pixel_vue4, left_knee_pixel_vue4, left_ankle_pixel_vue4)

%Functions

%function called to get the world coordinates of a specific joint in a
%specific frame
%This is used for Step 1
function joint_world_coords = get_joint_world_coords(frame, joint)
load('Subject4-Session3-Take4_mocapJoints.mat');
x = mocapJoints(frame,joint,1); %array of 12 X coordinates
y = mocapJoints(frame,joint,2); % Y coordinates
z = mocapJoints(frame,joint,3); % Z coordinates
conf = mocapJoints(frame,joint,4); %confidence values
joint_world_coords = [x,y,z,conf];
end

function vue2_joint_pixel_coords = get_vue2_joint_pixel_coords(frame, joint)
vue2_pmat = [-0.759251863679378,-0.649110013194741,0.0468273239649182,137.715365606944;-0.138093351133920,0.0903736400713044,-0.986287398050318,805.522696238255;0.635977070240741,-0.755307087184137,-0.158254131628272,7336.52052400539;0,0,0,1];
vue2_joint_pixel_coords = zeros(length(frame),2);%This preallocates the answer to increase the efficiency of the loop
%Iterate through the world coordinates for every frame
for i = 1:length(frame)
    %Get the world coordinates of joint in the frame.
    joint_world_coords = [joint(i,1); joint(i,2); joint(i,3); joint(i,4)];
    %Convert world coords to camera coords:
    joint_cam_coords = vue2_pmat * joint_world_coords;
    %Convert camera coords to film coords
    joint_film_coords = [1557.8 0 0 0; 0 1557.8 0 0; 0 0 1 0] * joint_cam_coords;
    joint_film_x = joint_film_coords(1) / joint_film_coords(3);
    joint_film_y = joint_film_coords(2) / joint_film_coords(3);
    %Convert film coords to pixel coords
    joint_pixel_x = joint_film_x + 976.0397;
    joint_pixel_y = joint_film_y + 562.8225;
    
    vue2_joint_pixel_coords(i,:) = [joint_pixel_x, joint_pixel_y];
end
end

function vue4_joint_pixel_coords = get_vue4_joint_pixel_coords(frame,joint)
vue4_pmat = [-0.816427158217721,0.576669940506233,0.0299745732425930,388.768652422159;0.103963925393142,0.197851965453977,-0.974703084012177,295.060869156745;-0.568012497698235,-0.792657794689825,-0.221484588574639,7283.00415138530; 0,0,0,1];
vue4_joint_pixel_coords = zeros(length(frame),2); %This preallocates the answer to increase the efficiency of the loop
%Iterate through the world coordinates for every frame
for i = 1:length(frame)
    %Get the world coordinates of joint in the frame.
    joint_world_coords = [joint(i,1); joint(i,2); joint(i,3); joint(i,4)];
    %Convert world coords to camera coords:
    joint_cam_coords = vue4_pmat * joint_world_coords;
    %Convert camera coords to film coords
    joint_film_coords = [1518.5 0 0 0; 0 1518.5 0 0; 0 0 1 0] * joint_cam_coords;
    joint_film_x = joint_film_coords(1) / joint_film_coords(3);
    joint_film_y = joint_film_coords(2) / joint_film_coords(3);
    %Convert film coords to pixel coords
    joint_pixel_x = joint_film_x + 975.0173;
    joint_pixel_y = joint_film_y + 554.8964;
    %Add the calculated points to the answer matrix
    vue4_joint_pixel_coords(i,:) = [joint_pixel_x, joint_pixel_y];
end
end

%This function calculates the pixel coordinates of the vue2 camera in vue4
function vue2_in_vue4 = get_v2_in_v4()
vue4_pmat = [-0.816427158217721,0.576669940506233,0.0299745732425930,388.768652422159;0.103963925393142,0.197851965453977,-0.974703084012177,295.060869156745;-0.568012497698235,-0.792657794689825,-0.221484588574639,7283.00415138530; 0,0,0,1];
vue2_world_coords = [-4450.06085208568;5557.92035162029;1949.06272680255;1];
%Convert world coords to camera coords:
vue2_cam_coords = vue4_pmat * vue2_world_coords;
%Convert camera coords to film coords
vue2_film_coords = [1518.5 0 0 0; 0 1518.5 0 0; 0 0 1 0] * vue2_cam_coords;
vue2_film_x = vue2_film_coords(1) / vue2_film_coords(3);
vue2_film_y = vue2_film_coords(2) / vue2_film_coords(3);
%Convert film coords to pixel coords
vue2_pixel_x = vue2_film_x + 975.0173;
vue2_pixel_y = vue2_film_y + 554.8964;
%Return the pixel coordinate
vue2_in_vue4 = [vue2_pixel_x, vue2_pixel_y];
end

function vue4_in_vue2 = get_v4_in_v2()
vue2_pmat = [-0.759251863679378,-0.649110013194741,0.0468273239649182,137.715365606944;-0.138093351133920,0.0903736400713044,-0.986287398050318,805.522696238255;0.635977070240741,-0.755307087184137,-0.158254131628272,7336.52052400539;0,0,0,1];
vue4_world_coords = [4.423562978688680e+03;5.490360440799710e+03;1.889016742748900e+03;1];
%Get the world coordinates of joint in the frame.
    %Convert world coords to camera coords:
    vue4_cam_coords = vue2_pmat * vue4_world_coords;
    %Convert camera coords to film coords
    vue4_film_coords = [1557.8 0 0 0; 0 1557.8 0 0; 0 0 1 0] * vue4_cam_coords;
    vue4_film_x = vue4_film_coords(1) / vue4_film_coords(3);
    vue4_film_y = vue4_film_coords(2) / vue4_film_coords(3);
    %Convert film coords to pixel coords
    vue4_pixel_x = vue4_film_x + 976.0397;
    vue4_pixel_y = vue4_film_y + 562.8225;
    %Return the pixel coordinate
    vue4_in_vue2 = [vue4_pixel_x, vue4_pixel_y];
end

%This function will get the vue2 camera viewing ray for a given joint 
function vue2_camera_ray = vue2_calc_viewing_ray(frame, joint_pixels)
%Calculate the vue2 u vector
vue2_camera_ray = vue2_transposed_rmat*inv(vue2.Kmat)*[joint_pixels(frame,1); joint_pixels(frame,2); 1];
end
%This function will get the vue4 camera viewing ray for a given joint 
function vue4_camera_ray = vue4_calc_viewing_ray(frame,joint_pixels)
%Calculate the vue4 u vector
vue4_camera_ray = vue4_transposed_rmat*inv(vue4.Kmat)*[joint_pixels(frame,1); joint_pixels(frame,2); 1];
end

function [a,b,d] = find_constants(u1,u2,u3,c1,c2)
%Calculate final c value
c_final = c2 - c1;
%Find constant values a,b,d
constants = inv([u1,u2,u3]) * c_final;
a = constants(1);
b = constants(2);
d = constants(3);
end

function recovered_3d_point = calc_triangulation(frame, joint_coords)
%vue2 camera parameters
vue2_Rmat = [-0.759251863679378,-0.649110013194741,0.0468273239649182;-0.138093351133920,0.0903736400713044,-0.986287398050318;0.635977070240741,-0.755307087184137,-0.158254131628272];
vue2_Kmat = [1.557796485665730e+03,0,9.760396728515630e+02;0,1.557796485665730e+03,5.628225097656250e+02;0,0,1];
vue2_transposed_rmat = transpose(vue2_Rmat);
vue2_t_vec = [137.715365606944;805.522696238255;7336.52052400539];
%vue4 camera parameters
vue4_Rmat = [-0.816427158217721,0.576669940506233,0.0299745732425930;0.103963925393142,0.197851965453977,-0.974703084012177;-0.568012497698235,-0.792657794689825,-0.221484588574639];
vue4_Kmat = [1.518513444384530e+03,0,9.750173339843750e+02;0,1.518513444384530e+03,5.548963623046880e+02;0,0,1];
vue4_transposed_rmat = transpose(vue4_Rmat);
vue4_t_vec = [388.768652422159;295.060869156745;7283.00415138530];
%Get c1 and c2 values, then get final c vector
vue2_c1 = -vue2_transposed_rmat*vue2_t_vec;
vue4_c2 = -vue4_transposed_rmat*vue4_t_vec;
%find u1, u2, and u3 vector values
u1 = vue2_transposed_rmat*inv(vue2_Kmat)*[joint_coords(frame,1); joint_coords(frame,2); 1];%vue2_calc_viewing_ray(frame, joint_coords);
u2 = vue4_transposed_rmat*inv(vue4_Kmat)*[joint_coords(frame,1); joint_coords(frame,2); 1];%vue4_calc_viewing_ray(frame, joint_coords);
u3_num = cross(u1,u2);
u3 = u3_num/norm(u3_num); 
%Find constant values a,b,d
[a,b,d] = find_constants(u1,u2,u3,vue2_c1,vue4_c2);
%Calculate p1 value
p1 = vue2_c1 + a*u1;
%Calculate p2 value
p2 = vue4_c2 + b*u2;
%Calculate final point value p
recovered_3d_point = p1+p2/2;
end

%These functions will display the points for each view for the input frame number
function verify_3Dto2D_vue2(frame,j1,j2,j3,j4,j5,j6,j7,j8,j9,j10,j11,j12)
mocapFnum = frame;
filenamevue2mp4 = 'Subject4-Session3-24form-Full-Take4-Vue2.mp4';
vue2video = VideoReader(filenamevue2mp4);

vue2video.CurrentTime = (mocapFnum-1)*(50/100)/vue2video.FrameRate;
vid2Frame = readFrame(vue2video);
image(vid2Frame)
axis on
hold on
j1n = j1(frame,:);
j2n = j2(frame,:);
j3n = j3(frame,:);
j4n = j4(frame,:);
j5n = j5(frame,:);
j6n = j6(frame,:);
j7n = j7(frame,:);
j8n = j8(frame,:);
j9n = j9(frame,:);
j10n = j10(frame,:);
j11n = j11(frame,:);
j12n = j12(frame,:);
plot(j1n(1),j1n(2), 'r.', 'Markersize', 10, 'Linewidth', 1);
plot(j2n(1),j2n(2), 'b.', 'Markersize', 10, 'Linewidth', 1);
plot(j3n(1),j3n(2), 'y.', 'Markersize', 10, 'Linewidth', 1);
plot(j4n(1),j4n(2), 'w.', 'Markersize', 10, 'Linewidth', 1);
plot(j5n(1),j5n(2), 'm.', 'Markersize', 10, 'Linewidth', 1);
plot(j6n(1),j6n(2), 'c.', 'Markersize', 10, 'Linewidth', 1);
plot(j7n(1),j7n(2), 'r.', 'Markersize', 10, 'Linewidth', 1);
plot(j8n(1),j8n(2), 'b.', 'Markersize', 10, 'Linewidth', 1);
plot(j9n(1),j9n(2), 'y.', 'Markersize', 10, 'Linewidth', 1);
plot(j10n(1),j10n(2), 'w.', 'Markersize', 10, 'Linewidth', 1);
plot(j11n(1),j11n(2), 'm.', 'Markersize', 10, 'Linewidth', 1);
plot(j12n(1),j12n(2), 'c.', 'Markersize', 10, 'Linewidth', 1);
%Plot the lines for the skeleton
plot([j1n(1),j4n(1)],[j1n(2),j4n(2)], 'r-', 'Markersize', 10, 'Linewidth', 1);
plot([j1n(1),j2n(1)],[j1n(2),j2n(2)], 'r-', 'Markersize', 10, 'Linewidth', 1);
plot([j2n(1),j3n(1)],[j2n(2),j3n(2)], 'r-', 'Markersize', 10, 'Linewidth', 1);
plot([j4n(1),j5n(1)],[j4n(2),j5n(2)], 'r-', 'Markersize', 10, 'Linewidth', 1);
plot([j5n(1),j6n(1)],[j5n(2),j6n(2)], 'r-', 'Markersize', 10, 'Linewidth', 1);
plot([j7n(1),j10n(1)],[j7n(2),j10n(2)], 'r-', 'Markersize', 10, 'Linewidth', 1);
plot([j7n(1),j8n(1)],[j7n(2),j8n(2)], 'r-', 'Markersize', 10, 'Linewidth', 1);
plot([j8n(1),j9n(1)],[j8n(2),j9n(2)], 'r-', 'Markersize', 10, 'Linewidth', 1);
plot([j10n(1),j11n(1)],[j10n(2),j11n(2)], 'r-', 'Markersize', 10, 'Linewidth', 1);
plot([j11n(1),j12n(1)],[j11n(2),j12n(2)], 'r-', 'Markersize', 10, 'Linewidth', 1);
spinexs = (j1n(1)+j4n(1))/2;
spineys = (j1n(2)+j4n(2))/2;
spinexh = (j7n(1)+j10n(1))/2;
spineyh = (j7n(2)+j10n(2))/2;
plot([spinexs,spinexh],[spineys,spineyh], 'r-', 'Markersize', 10, 'Linewidth', 1);
title(["Vue2 Frame", frame]);
end

function verify_3Dto2D_vue4(frame,j1,j2,j3,j4,j5,j6,j7,j8,j9,j10,j11,j12)
mocapFnum = frame;
filenamevue4mp4 = 'Subject4-Session3-24form-Full-Take4-Vue4.mp4';
vue4video = VideoReader(filenamevue4mp4);
vue4video.CurrentTime = (mocapFnum-1)*(50/100)/vue4video.FrameRate;
vid2Frame = readFrame(vue4video);
image(vid2Frame)
axis on
hold on
j1n = j1(frame,:);
j2n = j2(frame,:);
j3n = j3(frame,:);
j4n = j4(frame,:);
j5n = j5(frame,:);
j6n = j6(frame,:);
j7n = j7(frame,:);
j8n = j8(frame,:);
j9n = j9(frame,:);
j10n = j10(frame,:);
j11n = j11(frame,:);
j12n = j12(frame,:);
plot(j1n(1),j1n(2), 'r.', 'Markersize', 10, 'Linewidth', 1);
plot(j2n(1),j2n(2), 'b.', 'Markersize', 10, 'Linewidth', 1);
plot(j3n(1),j3n(2), 'y.', 'Markersize', 10, 'Linewidth', 1);
plot(j4n(1),j4n(2), 'w.', 'Markersize', 10, 'Linewidth', 1);
plot(j5n(1),j5n(2), 'm.', 'Markersize', 10, 'Linewidth', 1);
plot(j6n(1),j6n(2), 'c.', 'Markersize', 10, 'Linewidth', 1);
plot(j7n(1),j7n(2), 'r.', 'Markersize', 10, 'Linewidth', 1);
plot(j8n(1),j8n(2), 'b.', 'Markersize', 10, 'Linewidth', 1);
plot(j9n(1),j9n(2), 'y.', 'Markersize', 10, 'Linewidth', 1);
plot(j10n(1),j10n(2), 'w.', 'Markersize', 10, 'Linewidth', 1);
plot(j11n(1),j11n(2), 'm.', 'Markersize', 10, 'Linewidth', 1);
plot(j12n(1),j12n(2), 'c.', 'Markersize', 10, 'Linewidth', 1);
%Plot the lines for the skeleton
plot([j1n(1),j4n(1)],[j1n(2),j4n(2)], 'r-', 'Markersize', 10, 'Linewidth', 1);
plot([j1n(1),j2n(1)],[j1n(2),j2n(2)], 'r-', 'Markersize', 10, 'Linewidth', 1);
plot([j2n(1),j3n(1)],[j2n(2),j3n(2)], 'r-', 'Markersize', 10, 'Linewidth', 1);
plot([j4n(1),j5n(1)],[j4n(2),j5n(2)], 'r-', 'Markersize', 10, 'Linewidth', 1);
plot([j5n(1),j6n(1)],[j5n(2),j6n(2)], 'r-', 'Markersize', 10, 'Linewidth', 1);
plot([j7n(1),j10n(1)],[j7n(2),j10n(2)], 'r-', 'Markersize', 10, 'Linewidth', 1);
plot([j7n(1),j8n(1)],[j7n(2),j8n(2)], 'r-', 'Markersize', 10, 'Linewidth', 1);
plot([j8n(1),j9n(1)],[j8n(2),j9n(2)], 'r-', 'Markersize', 10, 'Linewidth', 1);
plot([j10n(1),j11n(1)],[j10n(2),j11n(2)], 'r-', 'Markersize', 10, 'Linewidth', 1);
plot([j11n(1),j12n(1)],[j11n(2),j12n(2)], 'r-', 'Markersize', 10, 'Linewidth', 1);
spinexs = (j1n(1)+j4n(1))/2;
spineys = (j1n(2)+j4n(2))/2;
spinexh = (j7n(1)+j10n(1))/2;
spineyh = (j7n(2)+j10n(2))/2;
plot([spinexs,spinexh],[spineys,spineyh], 'r-', 'Markersize', 10, 'Linewidth', 1);
title(["Vue4 Frame", frame]);
end