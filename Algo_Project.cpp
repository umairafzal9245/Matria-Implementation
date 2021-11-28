#include<iostream>
#include<algorithm>
#include<vector>
#include<fstream>
#include<iomanip>
#include<limits>
#include<string>
#include<list>
#include<queue>
#include <chrono>
#include<sstream>
using namespace std;
using namespace std::chrono;

struct Vertex
{
	int nodenumber;
	vector<pair<Vertex, double>> adjacencylist;
	Vertex()
	{
		nodenumber = 0;
	}
};

class Graph
{
public:
	int numberofvertices;
	vector<Vertex> vertexlist;
	pair<int, double>** matrix;
	Graph()
	{
		numberofvertices = 0;
		matrix = nullptr;
	}
	Graph& operator=(const Graph& obj)
	{
		numberofvertices = obj.numberofvertices;
		matrix = new pair<int, double> * [numberofvertices];
		for (int i = 0; i < numberofvertices; i++)
		{
			matrix[i] = new pair<int, double>[numberofvertices];
		}
		for (int i = 0; i < numberofvertices; i++)
		{
			for (int j = 0; j < numberofvertices; j++)
			{
				matrix[i][j].first = obj.matrix[i][j].first;
				matrix[i][j].second = obj.matrix[i][j].second;
			}
		}

		for (int i = 0; i < numberofvertices; i++)
		{
			Vertex* Node = new Vertex;
			Node->nodenumber = i;
			vertexlist.push_back(*Node);
		}

		for (int i = 0; i < numberofvertices; i++)
		{
			for (int j = 0; j < numberofvertices; j++)
			{
				if (matrix[i][j].second != 0)
				{
					Vertex Node;
					Node.nodenumber = j;
					vertexlist.at(i).adjacencylist.push_back(make_pair(Node, matrix[i][j].second));

				}
			}
		}
		return *this;
	}
	void AddVertex(int N)
	{
		numberofvertices = N;
		matrix = new pair<int, double> * [numberofvertices];
		for (int i = 0; i < numberofvertices; i++)
		{
			matrix[i] = new pair<int, double>[numberofvertices];
		}

		for (int i = 0; i < numberofvertices; i++)
		{
			for (int j = 0; j < numberofvertices; j++)
			{
				matrix[i][j].first = 0;
				matrix[i][j].second = 0;
			}
		}

		for (int i = 0; i < N; i++)
		{
			Vertex* Node = new Vertex;
			Node->nodenumber = i;
			vertexlist.push_back(*Node);
		}

		for (int i = 0; i < numberofvertices; i++)
		{
			for (int j = 0; j < numberofvertices; j++)
			{
				if (j == 0)
					matrix[i][j].first = i;
				else
					matrix[i][j].first = j;
			}
		}
	}
	double* gettokens(string line, int vertices)
	{
		double* array = new double[vertices];
		string* tokens = new string[vertices];
		int j = 0;
		for (int i = 0; line[i] != '\0'; i++)
		{

			if (line[i] == ',')
			{
				i++;
				j++;
			}
			else
				tokens[j] += line[i];
		}
		for (int i = 0; i < vertices; i++)
		{
			stringstream object(tokens[i]);
			object >> array[i];
			//array[i] = stod(tokens[i]);
		}
		return array;
	}
	void populategraph(string filename)
	{
		fstream object;
		object.open(filename, ios::in);
		string line;
		int vertices = 0;
		while (getline(object, line))
		{
			vertices++;
		}

		AddVertex(vertices);

		object.close();

		object.open(filename, ios::in);
		numberofvertices = vertices;

		for (int i = 0; i < vertices; i++)
		{
			getline(object, line);
			double* value = gettokens(line, vertices);
			for (int j = 0; j < vertices; j++)
			{
				matrix[i][j].second = value[j];
			}
		}

		for (int i = 0; i < vertices; i++)
		{
			for (int j = 0; j < vertices; j++)
			{
				if (matrix[i][j].second != 0)
				{
					Vertex Node;
					Node.nodenumber = j;
					vertexlist.at(i).adjacencylist.push_back(make_pair(Node, matrix[i][j].second));

				}
			}
		}
	}
	int Search(int node)
	{
		int i = 0;
		for (; i < vertexlist.size(); i++)
		{
			if (vertexlist.at(i).nodenumber == node)
				break;
		}
		return i;
	}
	vector<int> SearchMatrix(int node)
	{
		vector<int> indexes;
		for (int i = 0; i < numberofvertices; i++)
		{
			if (matrix[i][0].first == node)
			{
				indexes.push_back(i);
				break;
			}
		}
		for (int i = 0; i < numberofvertices; i++)
		{
			if (matrix[0][i].first == node)
			{
				indexes.push_back(i);
				break;
			}
		}

		return indexes;
	}
	void Display()
	{
		for (int i = 0; i < vertexlist.size(); i++)
		{
			cout << vertexlist.at(i).nodenumber;
			vector<pair<Vertex, double>> edgevector = vertexlist.at(i).adjacencylist;

			for (auto it = edgevector.begin(); it != edgevector.end(); it++)
			{
				cout << " -->";
				cout << " ( " << it->first.nodenumber << " , ";
				cout << fixed << setprecision(2) << it->second << " ) ";
			}
			cout << endl;
		}
		cout << endl;
		displaymatrix();
	}
	void displaymatrix()
	{
		for (int i = 0; i < numberofvertices; i++)
		{
			for (int j = 0; j < numberofvertices; j++)
			{
				cout << matrix[i][j].first << "," << fixed << setprecision(2) << matrix[i][j].second << "  ";
			}
			cout << endl;
		}
	}

	void removedependencies(int node)
	{
		int index = Search(node);// cout << "ind :"<< index<<endl;
		//if (index < 29) {
			vertexlist.at(index).adjacencylist.clear();
			vertexlist.erase(vertexlist.begin() + index);

			for (auto i = vertexlist.begin(); i != vertexlist.end(); i++)
			{
				vector<pair<Vertex, double>> innerlist = i->adjacencylist;
				int k = 0;
				for (auto j = innerlist.begin(); j != innerlist.end(); j++, k++)
				{
					if (j->first.nodenumber == node)
					{
						i->adjacencylist.erase(i->adjacencylist.begin() + k);
						break;
					}
				}
			}

			int row = SearchMatrix(node).at(0);
			int col = SearchMatrix(node).at(1);


			for (int i = 0; i < numberofvertices; i++)
			{
				matrix[row][i].second = 0;
				matrix[i][col].second = 0;
			}
		//}
	}

	double getweight(int src, int dest)
	{
		vector<int> srcpositions = SearchMatrix(src);
		vector<int> distpositions = SearchMatrix(dest);
		int row = srcpositions.at(0);
		int col = distpositions.at(1);
		return matrix[row][col].second;
	}
	bool bfs(int src, int last)
	{
		vector<bool> visited(numberofvertices, false);
		vector<int> vertices;
		vector<int> q;
		q.push_back(src);
		visited[src] = true;
		int vis;
		while (!q.empty())
		{
			vis = q[0];
			vertices.push_back(vis);
			q.erase(q.begin());
			for (int i = 0; i < numberofvertices; i++)
			{
				if (matrix[vis][i].second > 0 && (!visited[i]))
				{
					q.push_back(matrix[vis][i].first);
					visited[i] = true;
				}
			}
		}

		for (int i = 0; i < vertices.size(); i++)
		{
			if (vertices.at(i) == last)
				return true;
		}
		return false;
	}

	bool isConnected()
	{
		return bfs(0, numberofvertices - 1);
	}
	int minDistance(double dist[], bool sptSet[])
	{
		double min = DBL_MIN;
		int min_index;

		for (int v = 0; v < numberofvertices; v++)
			if (sptSet[v] == false && dist[v] >= min)
				min = dist[v], min_index = v;

		return min_index;
	}
	double moddijakstra(int sr)
	{
		int src = sr;
		double* dist = new double[numberofvertices];
		bool* sptSet = new bool[numberofvertices];

		for (int i = 0; i < numberofvertices; i++)
			dist[i] = DBL_MIN, sptSet[i] = false;


		dist[src] = 1;


		for (int count = 0; count < numberofvertices - 1; count++) {

			int u = minDistance(dist, sptSet);
			sptSet[u] = true;

			for (int v = 0; v < numberofvertices; v++)
				if (!sptSet[v] && matrix[u][v].second && dist[u] != DBL_MIN
					&& dist[u] * matrix[u][v].second > dist[v])
					dist[v] = dist[u] * matrix[u][v].second;
		}
		double sum = 0;
		for (int i = 0; i < numberofvertices; i++)
		{
			sum += dist[i];
		}
		return sum;
	}
	void printPath(int parent[], int j, string& path)
	{
		if (parent[j] == -1)
			return;

		printPath(parent, parent[j], path);
		path += to_string(j);
		//cout << "-->" << j;
	}
	string* printSolution(double dist[], int n, int parent[], int src, vector<int>& indexes)
	{
		string* path = new string[numberofvertices];
		path[src] = "0";
		for (int i = 0; i < indexes.size(); i++)
		{
			//cout << indexes.at(i) << "\t\t" << dist[indexes.at(i)] << "\t\t" << src;
			printPath(parent, indexes.at(i), path[indexes.at(i)]);
			//cout << endl;
		}
		for (int i = 0; i < indexes.size(); i++)
		{
			if (indexes.at(i) != src)
				path[indexes.at(i)].insert(0, 1, src);
		}
		return path;
	}
	int minDistancebet(double dist[], bool sptSet[], vector<int> indexes)
	{
		double min = DBL_MIN;
		int min_index = 0;

		for (int v = 0; v < indexes.size(); v++)
		{
			if (sptSet[indexes.at(v)] == false && dist[indexes.at(v)] >= min)
			{
				min = dist[indexes.at(v)], min_index = indexes.at(v);
			}
		}
		return min_index;
	}
	string* betweenesspathsusingdijakstra(int sr, vector<int> indexes)
	{
		int src = sr;
		double* dist = new double[numberofvertices];
		int* parent = new int[numberofvertices];
		bool* sptSet = new bool[numberofvertices];
		parent[src] = -1;

		for (int i = 0; i < indexes.size(); i++)
		{
			dist[indexes.at(i)] = DBL_MIN;
			sptSet[indexes.at(i)] = false;
		}
		dist[src] = 1;

		for (int count = 0; count < indexes.size() - 1; count++)
		{
			int u = minDistancebet(dist, sptSet, indexes);
			sptSet[u] = true;
			for (int v = 0; v < indexes.size(); v++)
			{
				if (!sptSet[indexes.at(v)] && dist[u] != DBL_MIN && matrix[u][indexes.at(v)].second && (dist[u] * matrix[u][indexes.at(v)].second > dist[indexes.at(v)]))
				{
					parent[indexes.at(v)] = u;
					dist[indexes.at(v)] = dist[u] * matrix[u][indexes.at(v)].second;
				}
			}
		}
		string* path = printSolution(dist, numberofvertices, parent, src, indexes);
		return path;
	}
	int BetweenessCentrality(int node)
	{
		vector<int> indexes;
		for (int i = 0; i < vertexlist.size(); i++)
		{
			indexes.push_back(vertexlist.at(i).nodenumber);
		}
		/*for (int i = 0; i < indexes.size(); i++)
		{
			cout << indexes.at(i) << " ";
		}
		cout << endl;
		*/
		string** pathmatrix = new string * [numberofvertices];
		for (int i = 0; i < numberofvertices; i++)
		{
			pathmatrix[i] = new string[numberofvertices];
		}

		for (int i = 0; i < indexes.size(); i++)
		{
			string* path = betweenesspathsusingdijakstra(indexes.at(i), indexes);
			for (int j = 0; j < indexes.size(); j++)
			{
				pathmatrix[indexes.at(i)][indexes.at(j)] += path[indexes.at(j)];
			}
		}

		string nod = to_string(node);
		int count = 0;
		for (int i = 0; i < indexes.size(); i++)
		{
			for (int j = 0; j < indexes.size(); j++)
			{
				if (indexes.at(i) != node && indexes.at(j) != node)
					if (pathmatrix[indexes.at(i)][indexes.at(j)].find(nod) != string::npos)
						count++;
			}
		}
		return count;
	}
	double DegreeCentrality(int node)
	{
		vector<int> indexes = SearchMatrix(node);
		double sum = 0;
		int row = indexes.at(0);
		for (int i = 0; i < numberofvertices; i++)
		{
			sum += matrix[row][i].second;
		}
		return sum;
	}
};


double ComputeCentralityCLO(Graph object, int node)
{
	return object.moddijakstra(node);
}

vector<int> ITERCENT_CLO(Graph object)
{
	Graph object2;
	object2 = object;
	vector<int> Set;
	do
	{
		double* values = new double[object2.numberofvertices];
		for (int i = 0; i < object2.numberofvertices; i++)
		{
			values[i] = ComputeCentralityCLO(object2, i);
		}
		double max = 0;
		for (int i = 0; i < object2.numberofvertices; i++)
		{
			if (values[i] > max)
				max = values[i];
		}
		int i = 0;
		for (; i < object2.numberofvertices; i++)
		{
			if (values[i] == max)
				break;
		}
		int removenode = i;
		Set.push_back(removenode);
		object2.removedependencies(removenode);
	} while (object2.isConnected());
	return Set;
}

double ComputeCentralityDEG(Graph object, int node)
{
	return object.DegreeCentrality(node);
}

vector<int> ITERCENT_DEG(Graph object)
{
	Graph object2;
	object2 = object;
	vector<int> Set;
	do
	{
		pair<int, double>* table = new pair<int, double>[object2.numberofvertices];
		for (int i = 0; i < object2.numberofvertices; i++)
		{
			table[i].first = i;
			table[i].second = ComputeCentralityDEG(object2, i);
			//cout << table[i].first << "  " << table[i].second << endl;
		}
		//cout << endl;
		pair<int, double> max;
		max.first = 0;
		max.second = 0;
		for (int i = 0; i < object2.numberofvertices; i++)
		{
			if (table[i].second > max.second)
			{
				max.second = table[i].second;
				max.first = table[i].first;
			}
		}
		//cout << "Removablenode = " << max.first << endl;
		Set.push_back(max.first);
		object2.removedependencies(max.first);
		//object2.Display();
	} while (object2.isConnected());
	/*cout << "set = ";
	for (int i = 0; i < Set.size(); i++)
	{
		cout << Set[i] << " ";
	}
	cout << endl;*/
	return Set;
}

double ComputeCentralityBET(Graph object, int node)
{
	return object.BetweenessCentrality(node);
}

vector<int> ITERCENT_BET(Graph object)
{
	Graph object2;
	object2 = object;
	vector<int> Set;

	do
	{
		pair<int, double>* table = new pair<int, double>[object2.numberofvertices];
		vector<int> indexes;
		for (int i = 0; i < object.vertexlist.size(); i++)
		{
			indexes.push_back(object.vertexlist.at(i).nodenumber);
		}
		for (int i = 0; i < indexes.size(); i++)
		{
			table[indexes.at(i)].first = indexes.at(i);
			table[indexes.at(i)].second = ComputeCentralityBET(object2, indexes.at(i));
			//cout << table[indexes.at(i)].first << "  " << table[indexes.at(i)].second << endl;
		}
		//cout << endl;
		pair<int, double> max;
		max.first = '0';
		max.second = 0;
		for (int i = 0; i < object2.numberofvertices; i++)
		{
			if (table[i].second > max.second)
			{
				max.second = table[i].second;
				max.first = table[i].first;
			}
		}
		//cout << "Removablenode = " << max.first << endl;
		Set.push_back(max.first);
		object2.removedependencies(max.first);
		//object2.Display();
	} while (object2.isConnected());
	/*cout << "set = ";
	for (int i = 0; i < Set.size(); i++)
	{
		cout << Set[i] << " ";
	}
	cout << endl;*/
	return Set;
}

void printvector(vector<int> object)
{
	for (int i = 0; i < object.size(); i++)
	{
		cout << object[i] << " ";
	}
	cout << endl;
}



vector<int> UNIFY(vector<int> SCLO, vector<int> SDEG, vector<int> SBET, Graph& obj)
{
	cout << "SCLO = ";
	printvector(SCLO);
	cout << "SDEG = ";
	printvector(SDEG);
	cout << "SBET = ";
	printvector(SBET);
	vector<int> finalresult;

	//UNIVERSAL AGREEMENTS
	vector<int> SCLOINTERSDEG;
	std::sort(SCLO.begin(), SCLO.end());
	std::sort(SDEG.begin(), SDEG.end());
	std::sort(SBET.begin(), SBET.end());

	std::set_intersection(SCLO.begin(), SCLO.end(),
		SDEG.begin(), SDEG.end(),
		back_inserter(SCLOINTERSDEG));
	std::set_intersection(SCLOINTERSDEG.begin(), SCLOINTERSDEG.end(),
		SBET.begin(), SBET.end(),
		back_inserter(finalresult));
	//UNIVERSAL AGREEMENTS END
	//SUPERNODE TRIADS

	vector<int> SCLO_U_SDEG;
	vector<int> SBET_MINUS_SCLOUSDEG;
	vector<int> SBET_U_SDEG;
	vector<int> SCLO_MINUS_SBETUSDEG;
	vector<int> SBET_U_SCLO;
	vector<int> SDEG_MINUS_SBETUSCLO;

	std::set_union(SCLO.begin(), SCLO.end(), SDEG.begin(), SDEG.end(), back_inserter(SCLO_U_SDEG));
	std::set_difference(SCLO_U_SDEG.begin(), SCLO_U_SDEG.end(), SBET.begin(), SBET.end(), back_inserter(SBET_MINUS_SCLOUSDEG));
	//printvector(SBET_MINUS_SCLOUSDEG);

	std::set_union(SBET.begin(), SBET.end(), SDEG.begin(), SDEG.end(), back_inserter(SBET_U_SDEG));
	std::set_difference(SBET_U_SDEG.begin(), SBET_U_SDEG.end(), SCLO.begin(), SCLO.end(), back_inserter(SCLO_MINUS_SBETUSDEG));
	//printvector(SCLO_MINUS_SBETUSDEG);

	std::set_union(SBET.begin(), SBET.end(), SCLO.begin(), SCLO.end(), back_inserter(SBET_U_SCLO));
	std::set_difference(SBET_U_SCLO.begin(), SBET_U_SCLO.end(), SDEG.begin(), SDEG.end(), back_inserter(SDEG_MINUS_SBETUSCLO));
	//printvector(SDEG_MINUS_SBETUSCLO);

	for (int x = 0; x < SBET_MINUS_SCLOUSDEG.size(); x++)
	{
		for (int y = 0; y < SCLO_MINUS_SBETUSDEG.size(); y++)
		{
			for (int z = 0; z < SDEG_MINUS_SBETUSCLO.size(); z++)
			{
				if (obj.getweight(SBET_MINUS_SCLOUSDEG[x], SCLO_MINUS_SBETUSDEG[y]) * obj.getweight(SCLO_MINUS_SBETUSDEG[y], SDEG_MINUS_SBETUSCLO[z]) * obj.getweight(SBET_MINUS_SCLOUSDEG[x], SDEG_MINUS_SBETUSCLO[z]) > 0)
				{
					if (!(std::find(finalresult.begin(), finalresult.end(), SBET_MINUS_SCLOUSDEG[x]) != finalresult.end()))
						finalresult.push_back(SBET_MINUS_SCLOUSDEG[x]);

					if (!(std::find(finalresult.begin(), finalresult.end(), SCLO_MINUS_SBETUSDEG[y]) != finalresult.end()))
						finalresult.push_back(SCLO_MINUS_SBETUSDEG[y]);

					if (!(std::find(finalresult.begin(), finalresult.end(), SDEG_MINUS_SBETUSCLO[z]) != finalresult.end()))
						finalresult.push_back(SDEG_MINUS_SBETUSCLO[z]);
				}
			}
		}
	}
	/*cout << "After adding supernode traids" << endl;
	printvector(finalresult);*/
	//SUPERNODE TRIADS END
	//SUPERNODE ADJACENCIES START

	vector<int> SCLO_INTER_SDEG;
	vector<int> SCLO_INTER_SDEG_MINUS_SBET;
	vector<int> SBET_INTER_SDEG;
	vector<int> SBET_INTER_SDEG_MINUS_SCLO;
	vector<int> SBET_INTER_SCLO;
	vector<int> SBET_INTER_SCLO_MINUS_SDEG;

	std::set_intersection(SCLO.begin(), SCLO.end(), SDEG.begin(), SDEG.end(), back_inserter(SCLO_INTER_SDEG));
	std::set_difference(SBET.begin(), SBET.end(), SCLO_INTER_SDEG.begin(), SCLO_INTER_SDEG.end(), back_inserter(SCLO_INTER_SDEG_MINUS_SBET));
	//printvector(SCLO_INTER_SDEG_MINUS_SBET);

	std::set_intersection(SBET.begin(), SBET.end(), SDEG.begin(), SDEG.end(), back_inserter(SBET_INTER_SDEG));
	std::set_difference(SCLO.begin(), SCLO.end(), SBET_INTER_SDEG.begin(), SBET_INTER_SDEG.end(), back_inserter(SBET_INTER_SDEG_MINUS_SCLO));
	//printvector(SBET_INTER_SDEG_MINUS_SCLO);

	std::set_intersection(SBET.begin(), SBET.end(), SCLO.begin(), SCLO.end(), back_inserter(SBET_INTER_SCLO));
	std::set_difference(SDEG.begin(), SDEG.end(), SBET_INTER_SCLO.begin(), SBET_INTER_SCLO.end(), back_inserter(SBET_INTER_SCLO_MINUS_SDEG));
	//printvector(SBET_INTER_SCLO_MINUS_SDEG);

	for (int x = 0; x < SBET.size(); x++)
	{
		for (int y = 0; y < SCLO_INTER_SDEG_MINUS_SBET.size(); y++)
		{
			if (obj.getweight(SBET[x], SCLO_INTER_SDEG_MINUS_SBET[y]) > 0)
			{
				if (!(std::find(finalresult.begin(), finalresult.end(), SBET[x]) != finalresult.end()))
				{
					finalresult.push_back(SBET[x]);
				}
				if (!(std::find(finalresult.begin(), finalresult.end(), SCLO_INTER_SDEG_MINUS_SBET[y]) != finalresult.end()))
				{
					finalresult.push_back(SCLO_INTER_SDEG_MINUS_SBET[y]);
				}
			}
		}
	}

	for (int x = 0; x < SCLO.size(); x++)
	{
		for (int y = 0; y < SBET_INTER_SDEG_MINUS_SCLO.size(); y++)
		{
			if (obj.getweight(SCLO[x], SBET_INTER_SDEG_MINUS_SCLO[y]) > 0)
			{
				if (!(std::find(finalresult.begin(), finalresult.end(), SCLO[x]) != finalresult.end()))
				{
					finalresult.push_back(SCLO[x]);
				}
				if (!(std::find(finalresult.begin(), finalresult.end(), SBET_INTER_SDEG_MINUS_SCLO[y]) != finalresult.end()))
				{
					finalresult.push_back(SBET_INTER_SDEG_MINUS_SCLO[y]);
				}
			}
		}
	}
	for (int x = 0; x < SDEG.size(); x++)
	{
		for (int y = 0; y < SBET_INTER_SCLO_MINUS_SDEG.size(); y++)
		{
			if (obj.getweight(SDEG[x], SBET_INTER_SCLO_MINUS_SDEG[y]) > 0)
			{
				if (!(std::find(finalresult.begin(), finalresult.end(), SDEG[x]) != finalresult.end()))
				{
					finalresult.push_back(SDEG[x]);
				}
				if (!(std::find(finalresult.begin(), finalresult.end(), SBET_INTER_SCLO_MINUS_SDEG[y]) != finalresult.end()))
				{
					finalresult.push_back(SBET_INTER_SCLO_MINUS_SDEG[y]);
				}
			}
		}
	}
	//cout << "After adding supernode adjacencies" << endl;
	//SUPERNODE ADJACENCIES END
	return finalresult;
}

void Matria(Graph object)
{
	vector<int> SCLO = ITERCENT_CLO(object);
	vector<int> SDEG = ITERCENT_DEG(object);
	vector<int> SBET = ITERCENT_BET(object);
	vector<int> Final = UNIFY(SCLO, SDEG, SBET, object);
	cout << "Final Result = "; printvector(Final);
}

int main()
{
	string filename = "";
	cout << "Enter the filename : ";
	cin >> filename;
	auto start = high_resolution_clock::now();
	Graph object, object2;
	object.populategraph(filename);
	Matria(object);
	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);
	cout << "Time taken by function: "
		<< duration.count() << " microseconds" << endl;
	return 0;
}