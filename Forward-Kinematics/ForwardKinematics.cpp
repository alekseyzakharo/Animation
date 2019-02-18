#include "ForwardKinematics.h"
#include <math.h>

#define PI 3.14159265359	

/*recursively calculate each joint's global position from the root (in-order traversal)*/
void ForwardKinematics::calculateJointPosWithQuaternion(Joint* joint)
{
	joint->GlobalPos = computeGlobalPosition(joint);

	Vector4 localQuat = computeLocalQuaternion(joint);

	joint->Globalquat = computeGlobalQuaternion(joint, localQuat);

	//recursively call self
	for (int i = 0; i < joint->childNum; i++)
	{
		calculateJointPosWithQuaternion(joint->getChild(i));
	}
}

//compute joint's global position with its local position and its parent joint's 
//global position and orientation
Vector4 ForwardKinematics::computeGlobalPosition(Joint* joint)
{
	Vector4 LocalGlobalPos = joint->LocalPos;
	if (joint->parent != NULL)
	{
		Vector4 invParentQuat(joint->parent->Globalquat.x*-1, joint->parent->Globalquat.y*-1, joint->parent->Globalquat.z*-1, joint->parent->Globalquat.w);
		LocalGlobalPos = quaternionMultiplication(quaternionMultiplication(joint->parent->Globalquat, Vector4(joint->LocalPos.x, joint->LocalPos.y, joint->LocalPos.z, 0)), invParentQuat);
		LocalGlobalPos.x = joint->parent->GlobalPos.x + LocalGlobalPos.x;
		LocalGlobalPos.y = joint->parent->GlobalPos.y + LocalGlobalPos.y;
		LocalGlobalPos.z = joint->parent->GlobalPos.z + LocalGlobalPos.z;
	}
	return LocalGlobalPos;
}


//get local rotation quaternion with frame data and rotation order (euler angle order)
Vector4 ForwardKinematics::computeLocalQuaternion(Joint* joint)
{
	float RadianX = joint->rot_angle1*PI / 180.0;
	float RadianY = joint->rot_angle2*PI / 180.0;
	float RadianZ = joint->rot_angle3*PI / 180.0;

	Vector4 vecX(sin(RadianX / 2), 0, 0, cos(RadianX / 2));
	Vector4 vecY(0, sin(RadianY / 2), 0, cos(RadianY / 2));
	Vector4 vecZ(0, 0, sin(RadianZ / 2), cos(RadianZ / 2));

	return quaternionMultiplication(quaternionMultiplication(quaternionMultiplication(Vector4(0, 0, 0, 1), vecX),vecY),vecZ);
}

//get global rotation quaternion with local rotation quaternion
Vector4 ForwardKinematics::computeGlobalQuaternion(Joint* joint, Vector4 localQuat)
{
	if (joint->parent != NULL)
		localQuat = quaternionMultiplication(joint->parent->Globalquat,localQuat);
	//default returnn
	return localQuat;
}

//calculate quaternion quat1*quat2
Vector4 ForwardKinematics::quaternionMultiplication(Vector4 quat1, Vector4 quat2)
{
	Vector3 v1 = Vector3(quat1.x, quat1.y, quat1.z);
	Vector3 v2 = Vector3(quat2.x, quat2.y, quat2.z);
	float x1 = quat1.x; float y1 = quat1.y; float z1 = quat1.z;
	float x2 = quat2.x; float y2 = quat2.y; float z2 = quat2.z;

	float w1 = quat1.w;
	float w2 = quat2.w;

	Vector3 v3 = w1 * v2 + w2 * v1 + v1.cross(v2);
	float w3 = w1 * w2 - v1.dot(v2);
	Vector4 result = Vector4(v3.x, v3.y, v3.z, w3);

	return result;
}