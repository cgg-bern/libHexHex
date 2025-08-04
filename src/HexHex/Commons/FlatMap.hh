#pragma once

#include <vector>
#include <utility>
#include <algorithm>
#include <stdexcept>

namespace HexHex
{

/**
 * @brief Generic Vector of Pairs that can be used like a map
 */
template <typename K, typename V>
class FlatMap
{
public:
    FlatMap() {}
    //~FlatMap() {}

    void set(const K& key, const V& value)
    {
        auto it = find(key);
        if (it != pairs.end()) {it->second = value;}
        else {pairs.emplace_back(key, value);}
    }

    bool get(const K& key, V& value) const
    {
        auto it = find(key);
        if (it != pairs.end()) {value = it->second; return true;}
        return false;
    }

    inline void reserve(char size)
    {
        pairs.reserve(size);
    }

    const V& at(const K& key) const
    {
        auto it = find(key);
        if (it != pairs.end()) {return it->second;}
        throw std::runtime_error("FlatMap: Key not found");
    }

    inline const std::pair<K, V>& getByIndex(const int& i) const
    {
        return pairs[i];
    }

    inline bool contains(const K& key) const
    {
        return find(key) != pairs.end();
    }

    inline size_t size() const
    {
        return pairs.size();
    }

    inline bool empty() const
    {
        return pairs.empty();
    }

    inline std::pair<K,V>& operator[](int index)
    {
        return pairs[index];
    }

    V& operator[](const K& key)
    {
        auto it = find(key);
        if (it != pairs.end()) {return it->second;}

        pairs.emplace_back(key, V());
        return pairs.back().second;
    }

    const V& operator[](const K& key) const
    {
        auto it = find(key);
        if (it != pairs.end()) {return it->second;}
        throw std::runtime_error("FlatMap: Key not found");
    }

    inline void clear()
    {
        pairs.clear();
    }

private:
    std::vector<std::pair<K, V>> pairs;

    typename std::vector<std::pair<K, V>>::iterator find(const K& key) {
        return std::find_if(pairs.begin(), pairs.end(), [&key](const std::pair<K, V>& pair) {
            return pair.first == key;
        });
    }

    typename std::vector<std::pair<K, V>>::const_iterator find(const K& key) const {
        return std::find_if(pairs.begin(), pairs.end(), [&key](const std::pair<K, V>& pair) {
            return pair.first == key;
        });
    }
};

}

