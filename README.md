# Fault-control

This is a repo for the course Fault-tolerant and Robust Control.

The repository is structured as follows:

    ├── .gitignore         <- ignore file for git, ensures that we only store files we care about
    ├── Dockerfile         <- Dockerfile which contains a ROS noetic installation and included BlueROV software
    ├── README.md          <- The top-level README for developers using this project.
    │
    ├── media/             <- Contains pretty pictures for the README.md
    |
    ├── ros_ws/            <- Main folder for the ROS packages used in the course
    │   └── src/ (*)            
    │       ├── uuv_simulator  <- Underwater simulator package
    │       ├── bluerov2       <- Main ROS package for the BlueROV, contains the controllers, transforms, etc.
    │       └── bluerov2_gazebo <- BlueROV Gazebo packages
    │
    ├── scripts/           <- General scripts
    │   ├── build_and_create_docker.sh  <- Builds the image and creates a container
    |   ├── get_course_updates.sh       <- Downloads the latest changes from course respository
    |   ├── save_my_changes.sh          <- Uploads changes to personal fork (optional)
    │   └── start_container.sh          <- Starts the container
    │
    └── Training_Sessions/ (*)   <- Folders of exercise material for each lecture, and the exercise manual in a README.md
        ├── TS1
        ├── TS2
        ├── ...
        └── TS8


![](./media/m1_about_this_mac.png)