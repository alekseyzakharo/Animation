#include "InverseKinematics.h"
#include <vector>

const double EPSI = 0.000001;

//IK workflow
void InverseKinematics::IK() {

	//check if the endeffector and the target are close enough
	Vector3 endPos = Vector3(endEffector->GlobalPos.x, endEffector->GlobalPos.y, endEffector->GlobalPos.z);
	Vector3 tarPos = Vector3(target->x, target->y, target->z);

	int i = 0;
	while ((endPos - tarPos).length() > threshold && i < iterNum)
	{
		if (mode == 0)
		{
			CCDMode();
		}
		else if (mode == 1)
		{
			JacobianMode();
			forwardKinematcis->forwardKinematicsComputation(base);
		}
		endPos = Vector3(endEffector->GlobalPos.x, endEffector->GlobalPos.y, endEffector->GlobalPos.z);
		i++;
	}
}

//CCD Mode IK for right arm
void InverseKinematics::CCDMode()
{
	Joint *current = endEffector;
	Vector4  quat1;
	while (current != base)
	{
		Vector3 Pt = Vector3(target->x, target->y, target->z); //target node
		Vector3 Pe = Vector3(endEffector->GlobalPos.x, endEffector->GlobalPos.y, endEffector->GlobalPos.z); //endeffector node
		Vector3 Pc = Vector3(current->parent->GlobalPos.x, current->parent->GlobalPos.y, current->parent->GlobalPos.z); //current node base

		Vector3 PePc = (Pe - Pc).normalize(); //endeffector - current
		Vector3 PtPc = (Pt - Pc).normalize(); //target - current

		float angle = PePc.dot(PtPc); //cos(o)

		if (abs(angle) < 0.9999f) //if 1.0 then no rotation
		{
			Vector3 axis = PePc.cross(PtPc);
			axis.normalize();

			float turnRadian = (float)acos(angle);

			quat1 = forwardKinematcis->buildQuaternionRotationWithRad(turnRadian, axis.x, axis.y, axis.z);
			current->parent->Globalquat = forwardKinematcis->quaternionMultiplication(quat1,current->parent->Globalquat);
			forwardKinematcis->forwardKinematicsComputation(base);
		}
		current = current->parent;
	}
}

//Entire Right Arm Jacobian IK
void InverseKinematics::JacobianMode()
{
	float step = 0.01f; //step size used to minimize rotation angle of each link in arm, to minimize error

	Joint *current = endEffector;
	Eigen::MatrixXf mat(3,0);

	std::vector<Joint*> nodes;
	std::vector<Vector3> rotationAxis;

	while (current != base)
	{
		Vector3 Pt = Vector3(target->x, target->y, target->z); //target node
		Vector3 Pe = Vector3(endEffector->GlobalPos.x, endEffector->GlobalPos.y, endEffector->GlobalPos.z); //endeffector node
		Vector3 Pc = Vector3(current->parent->GlobalPos.x, current->parent->GlobalPos.y, current->parent->GlobalPos.z); //current node base

		Vector3 PePc = (Pe - Pc).normalize(); //endeffector - current
		Vector3 PtPc = (Pt - Pc).normalize(); //target - current

		//calculate axis = ai
		Vector3 axis = PePc.cross(PtPc);

		//save rotation axis for later
		rotationAxis.push_back(axis);

		//jacobian column (for each degree of freedom) = (rotation axis of node) x (vector to endeffector) == (ai*r)
		Vector3 column = axis.cross(PePc);

		//add a column to the matrix
		mat.conservativeResize(mat.rows(),mat.cols()+1);
		mat.col(mat.cols()-1) = Eigen::Vector3f(column.x, column.y, column.z).transpose();

		//save current node for figuring out qaternion later
		nodes.push_back(current);
		current = current->parent;
	}

	//multiply transpose
	mat.transposeInPlace();
	Eigen::MatrixXf dist(3, 1);

	//vp = vector from endeffector to target
	dist.col(0) = Eigen::Vector3f(target->x - endEffector->GlobalPos.x, target->y - endEffector->GlobalPos.y, target->z - endEffector->GlobalPos.z);

	//multiply jacobian transpose by vp to get angles for each link in arm
	mat *= dist;

	//for each node figure out its quaternion
	for (int i = 0; i < nodes.size(); i++)
	{
		Joint *j = nodes.at(i);
		Vector3 rot = rotationAxis.at(i);

		//multiply the angle calcualte by the stepsize to minimize error
		Vector4 quat = forwardKinematcis->buildQuaternionRotationWithRad(mat(i, 0)*step, rot.x, rot.y, rot.z);
		j->parent->Globalquat = forwardKinematcis->quaternionMultiplication(quat, j->parent->Globalquat);
	}
	//Do forward comutation only once after you figure out each of the nodes quaternion
	forwardKinematcis->forwardKinematicsComputation(base);
}