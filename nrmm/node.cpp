#include "node.h"

#include "node.h"

namespace npbayes
{
	void Node::GetLeaves(Set &leaves)
	{
		Queue ndque;
		ndque.push(this);
		while (!ndque.empty()) {
			Node *nd = ndque.front();
			ndque.pop();
			if (nd->IsLeaf())
				leaves.insert(nd);
			else {
				ndque.push(nd->ch0);
				ndque.push(nd->ch1);
			}
		}
	}

	void Node::GetLeaves(Queue &leaves)
	{
		Queue ndque;
		ndque.push(this);
		while (!ndque.empty()) {
			Node *nd = ndque.front();
			ndque.pop();
			if (nd->IsLeaf())
				leaves.push(nd);
			else {
				ndque.push(nd->ch0);
				ndque.push(nd->ch1);
			}
		}
	}

	Node * Node::GetRoot(void)
	{
		Node *rt = this;
		while (rt->par) rt = rt->par;
		return rt;
	}

	void Node::Shuffle(const Set &ndset, Queue &ndque)
	{
		std::vector<Node*> temp(ndset.begin(), ndset.end());
		std::random_shuffle(temp.begin(), temp.end());
		foreach(it, temp) ndque.push(it);
	}

	void Delete(Node *nd)
	{
		Node::Queue ndque;
		ndque.push(nd);
		while (!ndque.empty()) {
			nd = ndque.front();
			ndque.pop();
			if (!nd->IsLeaf()) {
				ndque.push(nd->ch0);
				ndque.push(nd->ch1);
			}
			delete nd;
		}
	}

	Node * CopySub(Node *src)
	{
		Node *dst = new Node(*src);
		if (src->ss) dst->ss = src->ss->Copy();
		return dst;
	}

	Node * Copy(Node *src)
	{
		Node *dst = CopySub(src);
		std::queue<std::pair<Node*, Node*>> ndque;
		ndque.push(std::pair<Node*, Node*>(src, dst));
		while (!ndque.empty()) {
			Node *src_temp = ndque.front().first;
			Node *dst_temp = ndque.front().second;
			ndque.pop();
			if (!src_temp->IsLeaf()) {
				dst_temp->ch0 = CopySub(src_temp->ch0);
				dst_temp->ch1 = CopySub(src_temp->ch1);
				dst_temp->ch0->par = dst_temp;
				dst_temp->ch1->par = dst_temp;				
				ndque.push(std::pair<Node*, Node*>(src_temp->ch0, dst_temp->ch0));
				ndque.push(std::pair<Node*, Node*>(src_temp->ch1, dst_temp->ch1));
			}
		}
		return dst;
	}
}