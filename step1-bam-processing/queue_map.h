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

#ifndef MODALIB_ORDERED_QUEUE_MAP_H
#define MODALIB_QUEUE_MAP_H

#include <queue>
#include <unordered_map>

template <class K, class V>
class OrderedQueueMap {
private:
    bool auto_pop;
    size_t max_size;
    std::queue<K> keys;
    std::unordered_map<K, V> values;
public:
    OrderedQueueMap() = delete;
    explicit OrderedQueueMap(size_t max_size, bool auto_pop) : auto_pop(auto_pop), max_size(max_size), keys(), values() {
        values.reserve(max_size);
    }
    explicit OrderedQueueMap(size_t max_size) : OrderedQueueMap(max_size, false) {}
    ~OrderedQueueMap() = default;
private:
    inline bool contains_value(const K &key) const { return this->values.find(key) != this->values.end(); }
public:
    void add(const K& key, const V& value) {

        // Remove oldest value if size is exceeded
        if (this->auto_pop) {
            if ((this->values.size() - 1 >= this->max_size) || (this->keys.size() - 1 >= 2 * this->max_size)) {
                this->pop();
            }
        }

        // Add the new entry
        this->keys.push(key);
        this->values[key] = value;
    }
    std::pair<K, V> front() const {
        const K& key = this->keys.front();
        const V& value = this->values.at(key);
        return {key, value};
    }
    inline void remove(const K& key) { this->values.erase(key); }
    inline bool contains(const K& key) { return this->values.find(key) != this->values.end(); }
    inline V& get(const K& key) { return this->values.at(key); }
    inline size_t size() const { return this->values.size(); }
    inline bool size_exceeded() const { return this->values.size() > this->max_size; }
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

#endif //MODALIB_ORDERED_QUEUE_MAP_H
