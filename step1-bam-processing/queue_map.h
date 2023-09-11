//  queue_map.h
//
//  Copyright 2023 Antoine Passemiers <antoine.passemiers@gmail.com>
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 2 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program; if not, write to the Free Software
//  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
//  MA 02110-1301, USA.

#ifndef MODALIB_QUEUE_MAP_H
#define MODALIB_QUEUE_MAP_H

#include <queue>
#include <unordered_map>

template <class K, class V>
class QueueMap {
private:
    size_t max_size;
    std::queue<K> keys;
    std::unordered_map<K, V> values;
public:
    QueueMap() = delete;
    QueueMap(size_t max_size) : max_size(max_size), keys(), values() {}
    ~QueueMap() = default;
    void add(const K& key, const V& value) {
        if (this->keys.size() - 1 >= this->max_size) this->pop();
        this->keys.push(key);
        this->values[key] = value;
    }
    inline void remove(const K& key) { this->values.erase(key); }
    inline bool contains(const K& key) { return this->values.find(key) != this->values.end(); }
    inline V& get(const K& key) { return this->values.at(key); }
    inline size_t size() const { return this->values.size(); }
private:
    void pop() {
        if (this->keys.size() == 0) {
            this->values.clear();
        } else {
            const K &key = this->keys.front();
            this->values.erase(key);
            this->keys.pop();
        }
    }
};

#endif //MODALIB_QUEUE_MAP_H
