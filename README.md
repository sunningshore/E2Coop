# $E^2Coop$ 

Energy efficiency is critical to trajectory planning for UAV swarms in obstacle detection and avoidance (ODA). This chapter presents a new scheme $E^2Coop$ that tightly couples two conventional methods: Artificial Potential Field (APF) with Particle Swarm Planning (PSO) in trajectory planning for ODA. In $E^2Coop$, swarm members perform trajectory planning cooperatively to avoid collisions in an energy-efficient manner. $E^2Coop$ exploits the advantages of the Active Contour Model in image processing for trajectory planning. Each swarm member plans its trajectories on the contours of APF to save energy and avoid collisions with obstacles. Swarm members that fall within the safeguard distance of each other plan their trajectories on different contours to avoid collisions. 

This work was published as 
```
@inproceedings{huang2021e2coop,
  title={E2Coop: Energy Efficient and Cooperative Obstacle Detection and Avoidance for UAV Swarms},
  author={Huang, Shuangyao and Zhang, Haibo and Huang, Zhiyi},
  booktitle={Proceedings of the International Conference on Automated Planning and Scheduling},
  volume={31},
  pages={634--642},
  year={2021}
}
```
If you find this code helpful, please cite 
```
@misc{UAV-Swarm-e2coop, 
  author = {Huang, Shuangyao}, 
  title = {E2Coop: Energy Efficient and Cooperative Obstacle Detection and Avoidance for UAV Swarms}, 
  year = {2022},
  publisher = {GitHub},
  journal = {GitHub repository},
  howpublished = {\url{https://github.com/sunningshore/E2Coop}},
}
```
