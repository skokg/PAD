
// MIT License
//
// Copyright (c) 2019 SuYuxi
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.


#ifndef KdTree_H
#define KdTree_H

#include <vector>
#include <queue>
#include <memory>
#include <cmath>
#include <algorithm>
#include <iostream>
//#include <limits>

namespace kdtree {
	using namespace std;

	typedef float PointType;
	typedef vector<PointType> Point; //Presents Point type by vector<int> like (x, y, z)

	struct KdTreeNode {
		KdTreeNode(Point _val = Point(), std::shared_ptr<KdTreeNode> _leftNode = nullptr, std::shared_ptr<KdTreeNode> _rightNode = nullptr)
				  : val(_val)
				  , leftNode(_leftNode)
				  , rightNode(_rightNode)
				  {}
		Point val;
		std::shared_ptr<KdTreeNode> leftNode;
		std::shared_ptr<KdTreeNode> rightNode;
	};

	//leftNode store points that point[splitDim] < node[splitDim]
	//rightNode store points that point[splitDim] >= node[splitDim]
	class KdTree {
	public:
		KdTree() : dimension(0), root(nullptr) {
		}

		std::shared_ptr<KdTreeNode> getKdTreeRoot() {
			return root;
		}

		bool buildKdTree(vector<Point>& points) { //All the points must be of same dimension like (x, y ,z) or (x, y) or what ever you like
			if(points.empty() | points[0].empty()) return false;
			dimension = points[0].size();
			root = buildHelper(points, 0, points.size() - 1, 0);
			return root != nullptr ? true : false; //return false if cannot build kd Tree by input points;
		}

		//print all Kd Tree's nodes layer by layer using breadth first search
		void printKdTree() {
			queue<std::shared_ptr<KdTreeNode>> q;
			q.emplace(root);
			std::shared_ptr<KdTreeNode> node;
			cout << "Points are:" << endl;;
			while(!q.empty())
			{
				node = q.front();
				q.pop();
				if(node != nullptr)
				{
					cout << "(";
					for_each(node->val.begin(), node->val.end() - 1, [](PointType& num) { cout << num << ", ";});
					cout << *(node->val.end() - 1) << ")" << endl;
					q.emplace(node->leftNode);
					q.emplace(node->rightNode);
				}
				else
				{
					cout << "(nullptr)" << endl;;
				}
			}
			cout << "End." << endl;;
		}

		std::shared_ptr<KdTreeNode> findMin(std::shared_ptr<KdTreeNode> node, int32_t dim, int32_t depth) //find the node with the minimum value on dim dimension from depth
		{
			if(root == nullptr) return nullptr;
			std::shared_ptr<KdTreeNode> minimum = node;
			findMinHelper(node, minimum, dim, depth);
			return minimum;
		}

		bool addNode(Point point) {
			if(root == nullptr || point.size() != (unsigned long)dimension) return false;
			std::shared_ptr<KdTreeNode> node = root;
			int32_t depth = 0;
			int32_t curDim;
			while(true)
			{
				curDim = depth % dimension;
				if(point[curDim] < node->val[curDim])
				{
					if(node->leftNode == nullptr)
					{
						node->leftNode = std::make_shared<KdTreeNode>(point);
						return true;
					}
					node = node->leftNode;
				}
				else
				{
					if(node->rightNode == nullptr)
					{
						node->rightNode = std::make_shared<KdTreeNode>(point);
						return true;
					}
					node = node->rightNode;
				}
				depth += 1;
			}

			return false;
		}

		void deleteNode(Point point) {
			if(root == nullptr || point.size() != (unsigned long)dimension) return;
			if(deleteNodeHelper(root, point, 0)) { root = nullptr; }
		}

		std::shared_ptr<KdTreeNode> getNode(Point point) {
			if(root == nullptr || point.size() != (unsigned long)dimension) return nullptr;
			std::shared_ptr<KdTreeNode> node = root;
			int32_t depth = 0;
			int32_t curDim;
			while(node != nullptr)
			{
				if(node->val == point) return node;
				curDim = depth % dimension;
				if(node->val[curDim] > point[curDim])
				{
					node = node->leftNode;
				}
				else
				{
					node = node->rightNode;
				}
				depth += 1;
			}

			return nullptr;
		}

		std::shared_ptr<KdTreeNode> findNearestNode(Point point) {
			if(root == nullptr || point.size() != (unsigned long)dimension) return nullptr;
			std::shared_ptr<KdTreeNode> nearestNode = root;
			float minDist = calDist(point, nearestNode->val);
			findNearestNodeHelper(root, point, minDist, nearestNode, 0);
			return nearestNode;
		}

		vector<std::shared_ptr<KdTreeNode>> findNearestNodeCluster(Point point, float distance) { //find all nodes with the distance from which to the "point" is less than or equal to the "distance."
			vector<std::shared_ptr<KdTreeNode>> cluster;
			if(root == nullptr || point.size() != (unsigned long)dimension) return cluster;
			findNearestNodeClusterHelper(root, point, distance, cluster, 0);
			return cluster;
		}

	private:
		std::shared_ptr<KdTreeNode> buildHelper(vector<Point>& points, int32_t leftBorder, int32_t rightBorder, int32_t depth) { //recursively create kd Tree
			if(leftBorder > rightBorder) return nullptr;
			int32_t curDim = depth % dimension;
			int32_t midInx = leftBorder + (rightBorder - leftBorder) / 2;
			//sort(points.begin() + leftBorder, points.begin() + rightBorder + 1, [curDim](const Point& a, const Point& b) { return a[curDim] < b[curDim]; });

			// ----------------------
			// ni nujno uporabiti celotnega sorta - namesto tega se lahko uporabi nth_element z še nekaj dodatnega dela - to je potem kar hitrejse - upam da zadeva dela prav
			nth_element(points.begin() + leftBorder, points.begin() + midInx,  points.begin() + rightBorder + 1, [curDim](const Point& a, const Point& b) { return a[curDim] < b[curDim]; });
			// sedaj je treba še na levi strani skupaj prestaviti vrednosti, ki so enake kot midInx (nth_element funkcija tega ne garantira)
			if (midInx - leftBorder > 1)
				{
				int32_t last_taken_indx=midInx;
				if (points[midInx-1][curDim] == points[midInx][curDim])
					last_taken_indx--;

				for (long il=midInx-2; il >= leftBorder; il--)
					if (points[il][curDim] == points[midInx][curDim])
						{
						swap(points[il],points[last_taken_indx-1]);
						last_taken_indx--;
						}
				}
			// ----------------------



			while(midInx > leftBorder && points[midInx - 1][curDim] == points[midInx][curDim]) // keep all the points with point[splitDim] >= midInx[splitDim] on the right of midInx
			{
				midInx -= 1;
			}
			std::shared_ptr<KdTreeNode> node = std::make_shared<KdTreeNode>(points[midInx]);
			node->leftNode = buildHelper(points, leftBorder, midInx - 1, depth + 1);
			node->rightNode = buildHelper(points, midInx + 1, rightBorder, depth + 1);
			return node;
		}

		void findMinHelper(std::shared_ptr<KdTreeNode>& node, std::shared_ptr<KdTreeNode>& minimum, const int32_t& dim, const int32_t& depth) {
			if(node == nullptr) return;
			if(node->val[dim] < minimum->val[dim]) { minimum = node; }
			findMinHelper(node->leftNode, minimum, dim, depth + 1);
			if(depth % dimension != dim) { findMinHelper(node->rightNode, minimum, dim, depth + 1); }
		}

		bool deleteNodeHelper(std::shared_ptr<KdTreeNode>& node, const Point& point, const int32_t& depth) { //return true means the child node should be deleted
			if(node == nullptr) return false;
			int32_t curDim = depth % dimension;
			if(point == node->val)
			{
				if(node->rightNode != nullptr)
				{
					std::shared_ptr<KdTreeNode> minimumNode = findMin(node->rightNode, curDim, depth + 1);
					node->val = minimumNode->val; // do not swap(node->val, minimumNode) which would break the structure of kd-tree and result in not finding the node to delete
					if(deleteNodeHelper(node->rightNode, minimumNode->val, depth + 1)) { node->rightNode = nullptr; }
				}
				else if(node->leftNode != nullptr)
				{
					std::shared_ptr<KdTreeNode> minimumNode = findMin(node->leftNode, curDim, depth + 1);
					node->val = minimumNode->val;
					if(deleteNodeHelper(node->leftNode, minimumNode->val, depth + 1)) { node->leftNode = nullptr; }
					node->rightNode = node->leftNode;
					node->leftNode = nullptr;
				}
				else
				{
					return true; //return true to inform outter deleteNodeHelper to release this pointer;
				}
			}
			else
			{
				if(point[curDim] < node->val[curDim])
				{
					if(deleteNodeHelper(node->leftNode, point, depth + 1)) { node->leftNode = nullptr; }
				}
				else
				{
					if(deleteNodeHelper(node->rightNode, point, depth + 1)) { node->rightNode = nullptr; }
				}
			}
			return false;
		}

		void findNearestNodeHelper(std::shared_ptr<KdTreeNode>& node, const Point& point, float& minDist, std::shared_ptr<KdTreeNode>& nearestNode, const int32_t& depth) {
			if(node == nullptr) return;
			int32_t curDim = depth % dimension;
			float dist = calDist(point, node->val);
			if(dist < minDist)
			{
				nearestNode = node;
				minDist = dist;
			}
			if(point[curDim] < node->val[curDim])
			{
				findNearestNodeHelper(node->leftNode, point, minDist, nearestNode, depth + 1);
				if(node->val[curDim] - point[curDim] <= minDist)
				{
					findNearestNodeHelper(node->rightNode, point, minDist, nearestNode, depth + 1);
				}
			}
			else
			{
				findNearestNodeHelper(node->rightNode, point, minDist, nearestNode, depth + 1);
				if(point[curDim] - node->val[curDim] < minDist)
				{
					findNearestNodeHelper(node->leftNode, point, minDist, nearestNode, depth + 1);
				}
			}
		}

		void findNearestNodeClusterHelper(std::shared_ptr<KdTreeNode>& node, const Point& point, const float& distance, vector<std::shared_ptr<KdTreeNode>>& cluster, const int32_t& depth) {
			if(node == nullptr) return;
			int32_t curDim = depth % dimension;
			float dist = calDist(point, node->val);
			if(dist < distance)
			{
				cluster.emplace_back(node);
			}

			if(point[curDim] < node->val[curDim])
			{
				findNearestNodeClusterHelper(node->leftNode, point, distance, cluster, depth + 1);
				if(node->val[curDim] - point[curDim] <= distance)
				{
					findNearestNodeClusterHelper(node->rightNode, point, distance, cluster, depth + 1);
				}
			}
			else
			{
				findNearestNodeClusterHelper(node->rightNode, point, distance, cluster, depth + 1);
				if(point[curDim] - node->val[curDim] < distance)
				{
				findNearestNodeClusterHelper(node->leftNode, point, distance, cluster, depth + 1);
				}
			}
		}

		/*PointType calDist(const Point& a, const Point& b) {
			PointType sum = 0;
			for(unsigned int32_t inx = 0; inx < min(a.size(), b.size()); inx += 1)
			{
				sum += pow(b[inx] - a[inx], 2);
			}
			return sqrtf(sum);
		}*/

		float calDist(const Point& a, const Point& b) {
			PointType sum = 0;
			for(uint32_t inx = 0; inx < min(a.size(), b.size()); inx += 1)
			{
				float temp = (float)(b[inx] - a[inx]);
				sum += temp*temp;
			}
			return sqrtf(sum);
		}

		int32_t dimension;
		std::shared_ptr<KdTreeNode> root;

	};
}
#endif
